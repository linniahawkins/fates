program MLPhenology

  use               FatesConstantsMod, only : r8 => fates_r8
  use              FatesArgumentUtils, only : command_line_arg
  use             FatesUnitTestIOMod,  only : OpenNCFile, GetVar, CloseNCFile, RegisterNCDims

  implicit none

  ! define ML phenoogy pytorch model
  character(len=256) :: the_torch_model = "/glade/u/home/linnia/MLphenology/models/example_LSTM_model_lh.pt"
  character(len=256) :: the_tft_torch_model = "/glade/u/home/ayal/phenology-ml-clm/models/tft_scripted.pt"
  !character(len=256) :: the_torch_model = "/glade/u/home/ayal/phenology-ml-clm/models/example_LSTM_model_v1.pt"
  
  real(8), dimension(10) :: dummy_lai
  integer :: sos_flag, n

  character(len=:),                  allocatable :: datm_file            ! input DATM 
  real(r8),                          allocatable :: ta_min(:)         ! daily min air temperature [degC]
  real(r8),                          allocatable :: ta_max(:)         ! daily max air temperature [degC]
  real(r8),                          allocatable :: pr(:)            ! daily precipitation [mm]
  real(r8),                          allocatable :: sw(:)                ! daily shortwave radiation (W/m2)
  real(r8),                          allocatable :: lai(:)              ! daily LAI (m2/m2)
  real(r8),                          allocatable :: soilm(:)              ! daily soil moisture at layer 3 (kg/m2)
  real(r8),                          allocatable :: doy(:)              ! day of year
  real(r8),                          allocatable :: photo(:)              ! daily photoperiod (seconds)

  real(r8)                                       :: doy_arr(10)              ! DOY array
  real(r8)                                       :: out_data(1,5)       ! output from the lstm model (lai)
  real(r8)                                       :: out_data_tft(1,10)       ! output from the tft model (lai)
  real(r8)                                       :: dayofyear ! day of year 
  real(r8)                                       :: soilt            ! soil temperature at 12cm 
  real(r8)                                       :: onset_gdd      ! onset growing degree days 
  real(r8)                                       :: onset_gddflag  ! Onset freeze flag

  logical                                        :: do_onset       ! Flag if onset should happen

  ! Constants
  integer,  parameter :: N_HIST = 60          ! Days of historical LAI
  integer,  parameter :: N_PRED = 10          ! Days of predicted LAI  
  integer,  parameter :: N_TOTAL = 70         ! Total window size
  integer,  parameter :: TODAY_IDX = 60       ! Index of "today" in window
  integer,  parameter :: TOLERANCE = 5        ! ±days around today for detection

  real(r8), parameter :: SOS_DERIV_THRESH = 0.03_r8   ! Min derivative for SOS (LAI/day)
  real(r8), parameter :: EOS_DERIV_THRESH = -0.03_r8  ! Max derivative for EOS (LAI/day)

  ! Phenological states
  integer, parameter :: STATE_DORMANT = 0
  integer, parameter :: STATE_GROWING = 1

  ! Seasonal boundaries (day of year)
  integer, parameter :: WINTER_SOLSTICE = 355   ! ~Dec 21 (or 1 after new year)
  integer, parameter :: SUMMER_SOLSTICE = 172   ! ~June 21
  
  ! Variables for phenology loop
  integer  :: n_days                          ! Total number of days in forcing
  integer  :: day_idx                         ! Current day index in loop
  integer  :: current_doy                     ! Current day of year
  integer  :: current_state                   ! Current phenological state
  integer  :: last_doy                        ! Previous day of year (for year detection)
  integer  :: year_idx                        ! Current year index
  real(r8) :: lai_70day(N_TOTAL)              ! Combined 70-day LAI window
  logical  :: trigger_sos, trigger_eos        ! Transition flags

  ! Output arrays for detected phenology dates
  integer, parameter :: MAX_YEARS = 20
  integer  :: sos_doy_detected(MAX_YEARS)     ! Detected SOS day of year
  integer  :: eos_doy_detected(MAX_YEARS)     ! Detected EOS day of year
  
  ! File I/O
  integer  :: out_unit, lai_unit
  character(len=256) :: output_dir, phenology_file, lai_file

  ! Load forcing data
  datm_file = command_line_arg(1) ! one year of daily ta, pr, sw, lai
  call load_met_forcing(datm_file, ta_min, ta_max, pr, sw, lai, soilm, doy, photo)
  n_days = size(lai)

  print *, "========================================"
  print *, "ML Phenology Detection Test"
  print *, "========================================"
  print *, "Loaded forcing data with ", n_days, " days"
  print *, "LAI range: ", minval(lai), " to ", maxval(lai)
  print *, ""
  ! ======================================
  ! Test 1: Original CLM SeasonalDecidOnset function (for comparison)
  ! ======================================
  print *, "--- Test 1: CLM GDD-based onset (single day test) ---"
  dayofyear = 1.0_r8
  onset_gdd = 0.0_r8
  onset_gddflag = 1.0_r8
  soilt = ta_min(1) - 10.0_r8

  do_onset = SeasonalDecidOnset( onset_gdd, onset_gddflag, soilt, dayofyear )
  print *, "  Day 1: onset_gdd =", onset_gdd, ", onset_gddflag =", onset_gddflag
  print *, ""

! ======================================
  ! Test 2: ML-based phenology detection loop
  ! ======================================
  print *, "--- Test 2: ML-based phenology detection (using observed LAI as 'prediction') ---"
  print *, "Looping through forcing data to detect SOS and EOS..."
  print *, ""

  ! Initialize
  current_state = STATE_DORMANT
  year_idx = 1
  last_doy = 0
  sos_doy_detected = 0
  eos_doy_detected = 0

  ! Loop through days, starting at day 61 (need 60 days of history)
  ! and ending 10 days before the end (need 10 days for "prediction")
  do day_idx = N_HIST + 1, n_days - N_PRED

    current_doy = int(doy(day_idx))

    ! Detect year transition (DOY goes from 365/366 back to 1)
    if (current_doy < last_doy) then
      year_idx = year_idx + 1
      print *, "  --- Year", year_idx, "started ---"
      ! Reset state at start of new year if still growing (safety reset)
      if (current_state == STATE_GROWING) then
        print *, "  WARNING: Still in GROWING state at year start, resetting to DORMANT"
        current_state = STATE_DORMANT
      end if
    end if
    last_doy = current_doy

    ! Check array bounds
    if (year_idx > MAX_YEARS) exit

    ! ============================================
    ! Build the 70-day LAI window
    ! ============================================
    ! Historical: 60 days of observed LAI (days day_idx-60 to day_idx-1)
    lai_70day(1:N_HIST) = lai(day_idx - N_HIST : day_idx - 1)

    ! "Predicted": For now, use observed LAI as a stand-in for TFT predictions
    ! In full implementation, this would be: call run_tft_model(..., out_data_tft)
    ! and then: lai_70day(N_HIST+1:N_TOTAL) = out_data_tft(1,:)
    lai_70day(N_HIST+1:N_TOTAL) = lai(day_idx : day_idx + N_PRED - 1)

    ! ============================================
    ! Detect phenology transitions
    ! ============================================
    call detect_phenology_transition(lai_70day, current_doy, current_state, &
                                      trigger_sos, trigger_eos)

    ! Handle SOS trigger
    if (trigger_sos) then
      sos_doy_detected(year_idx) = current_doy
      current_state = STATE_GROWING
      print *, "  ** SOS triggered on DOY", current_doy, "(year", year_idx, ")"
      print *, "     LAI at trigger:", lai(day_idx)
    end if

    ! Handle EOS trigger
    if (trigger_eos) then
      eos_doy_detected(year_idx) = current_doy
      current_state = STATE_DORMANT
      print *, "  ** EOS triggered on DOY", current_doy, "(year", year_idx, ")"
      print *, "     LAI at trigger:", lai(day_idx)
    end if

  end do

  ! ======================================
  ! Summary of detected phenology
  ! ======================================
  print *, ""
  print *, "========================================"
  print *, "Summary of Detected Phenology Events"
  print *, "========================================"
  do n = 1, year_idx
    if (sos_doy_detected(n) > 0 .or. eos_doy_detected(n) > 0) then
      print *, "Year", n, ": SOS =", sos_doy_detected(n), ", EOS =", eos_doy_detected(n)
    end if
  end do
  print *, ""

  ! ======================================
  ! Write results to files
  ! ======================================
  output_dir = "ml_phenology_output"
  call system('mkdir -p ' // trim(output_dir))
  
  ! Write phenology events
  phenology_file = trim(output_dir) // '/phenology_events.csv'
  open(newunit=out_unit, file=trim(phenology_file), status='replace', action='write')
  write(out_unit, '(A)') 'year,sos_doy,eos_doy'
  do n = 1, year_idx
    write(out_unit, '(I0,",",I0,",",I0)') n, sos_doy_detected(n), eos_doy_detected(n)
  end do
  close(out_unit)
  print *, "Phenology events written to: ", trim(phenology_file)
  
  ! Write LAI time series
  lai_file = trim(output_dir) // '/lai_timeseries.csv'
  open(newunit=lai_unit, file=trim(lai_file), status='replace', action='write')
  write(lai_unit, '(A)') 'day,doy,lai'
  do n = 1, n_days
    write(lai_unit, '(I0,",",F10.3,",",F12.6)') n, doy(n), lai(n)
  end do
  close(lai_unit)
  print *, "LAI time series written to: ", trim(lai_file)
  print *, ""


  ! ====================================
  ! test get_SOS
  n = 10
  dummy_lai = (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0/)
  call get_sos(dummy_lai, n, sos_flag)
  print *, "start of season Flag:", sos_flag


  ! ========================================
  ! test tft

  call next_ten_days(dayofyear, doy_arr)
  call run_tft_model(the_tft_torch_model, ta_min, pr, sw, lai, doy_arr, out_data_tft)
  print *, "TFT predicted LAI now:", out_data_tft

  contains
  
    !-----------------------------------------------------------------------
    subroutine load_met_forcing ( datm_file, ta_min, ta_max, pr, sw, lai, soilm, doy, photo)
      ! 
      use FatesConstantsMod, only: r8 => fates_r8
      use FatesUnitTestIOMod, only: OpenNCFile, GetVar, CloseNCFile
    
      implicit none
    
      ! Arguments
      character(len=*), intent(in) :: datm_file
      real(r8), allocatable, intent(out) :: ta_min(:), ta_max(:), pr(:), sw(:), lai(:), soilm(:), doy(:), photo(:)
    
      ! Local
      integer :: ncid
    
      ! Allocate arrays
      allocate(ta_min(5844), ta_max(5844), pr(5844), sw(5844), lai(5844), soilm(5844), doy(5844), photo(5844))
    
      ! Open and read
      call OpenNCFile(trim(datm_file), ncid, 'read')
      
      call GetVar(ncid, 'ta_min', ta_min)
      call GetVar(ncid, 'ta_max', ta_max)
      call GetVar(ncid, 'pr', pr)
      call GetVar(ncid, 'sw', sw)
      call GetVar(ncid, 'lai', lai)
      call GetVar(ncid, 'soilm', soilm)
      call GetVar(ncid, 'doy', doy)
      call GetVar(ncid, 'photo', photo)
    
      call CloseNCFile(ncid)
    
    end subroutine load_met_forcing
    
    ! ----------------------------------------------------------------------
    subroutine compute_first_derivative(lai, n, deriv)
      !
      ! Compute first derivative of LAI using central finite differences.
      ! Uses forward difference at start and backward difference at end.
      !
      ! Input:
      !   lai(n)  - LAI time series
      !   n       - Length of array
      !
      ! Output:
      !   deriv(n) - First derivative (dLAI/dt, units: LAI/day)
      !
      !-----------------------------------------------------------------------
      
      implicit none
      
      ! Arguments
      integer,  intent(in)  :: n
      real(r8), intent(in)  :: lai(n)
      real(r8), intent(out) :: deriv(n)
      
      ! Local variables
      integer :: i
      
      ! Forward difference for first point
      deriv(1) = lai(2) - lai(1)
      
      ! Central differences for interior points
      do i = 2, n-1
        deriv(i) = (lai(i+1) - lai(i-1)) / 2.0_r8
      end do
      
      ! Backward difference for last point
      deriv(n) = lai(n) - lai(n-1)
      
    end subroutine compute_first_derivative
    !-----------------------------------------------------------------------
    subroutine find_local_extremum(deriv, n, search_start, search_end, &
                                    find_max, extremum_idx, extremum_val, found)
      !
      ! Find local maximum or minimum within a search window.
      !
      ! Input:
      !   deriv(n)     - First derivative array
      !   n            - Length of array
      !   search_start - Start index of search window
      !   search_end   - End index of search window
      !   find_max     - If .true., find local max; if .false., find local min
      !
      ! Output:
      !   extremum_idx - Index of the extremum (0 if not found)
      !   extremum_val - Value at the extremum
      !   found        - .true. if a valid local extremum was found
      !
      !-----------------------------------------------------------------------
      
      implicit none
      
      ! Arguments
      integer,  intent(in)  :: n
      real(r8), intent(in)  :: deriv(n)
      integer,  intent(in)  :: search_start
      integer,  intent(in)  :: search_end
      logical,  intent(in)  :: find_max
      integer,  intent(out) :: extremum_idx
      real(r8), intent(out) :: extremum_val
      logical,  intent(out) :: found
      
      ! Local variables
      integer  :: i
      integer  :: s_start, s_end
      logical  :: is_extremum
      
      ! Initialize outputs
      found = .false.
      extremum_idx = 0
      extremum_val = 0.0_r8
      
      ! Ensure search bounds are valid
      s_start = max(2, search_start)        ! Need neighbors for local check
      s_end   = min(n-1, search_end)
      
      if (s_start > s_end) return
      
      ! Search for local extremum
      do i = s_start, s_end
        
        if (find_max) then
          ! Check for local maximum: deriv(i) >= deriv(i-1) AND deriv(i) >= deriv(i+1)
          ! Use >= to handle flat peaks where consecutive values are equal
          is_extremum = (deriv(i) >= deriv(i-1)) .and. (deriv(i) >= deriv(i+1))
        else
          ! Check for local minimum: deriv(i) <= deriv(i-1) AND deriv(i) <= deriv(i+1)
          ! Use <= to handle flat valleys where consecutive values are equal
          is_extremum = (deriv(i) <= deriv(i-1)) .and. (deriv(i) <= deriv(i+1))
        end if
        
        if (is_extremum) then
          ! If this is the first extremum, or a more extreme value
          if (.not. found) then
            found = .true.
            extremum_idx = i
            extremum_val = deriv(i)
          else
            ! Keep the most extreme value in the window
            if (find_max .and. deriv(i) > extremum_val) then
              extremum_idx = i
              extremum_val = deriv(i)
            else if (.not. find_max .and. deriv(i) < extremum_val) then
              extremum_idx = i
              extremum_val = deriv(i)
            end if
          end if
        end if
        
      end do
    
    end subroutine find_local_extremum

    !-----------------------------------------------------------------------
    subroutine detect_phenology_transition(lai_70day, current_doy, current_state, &
                                            trigger_sos, trigger_eos)
      !
      ! Detect phenological transitions (SOS or EOS) based on the 70-day LAI window.
    !
    ! Method:
    !   1. Compute first derivative across the 70-day window
    !   2. Search for local max (SOS) or min (EOS) within ±5 days of "today" (day 60)
    !   3. Apply seasonal constraints:
    !      - SOS only allowed DOY 1-171 (winter solstice to summer solstice)
    !      - EOS only allowed DOY 172-365 (summer solstice to winter solstice)
    !   4. Apply state constraints:
    !      - SOS only if currently DORMANT
    !      - EOS only if currently GROWING
    !   5. Apply magnitude thresholds
    !
    ! Input:
    !   lai_70day(70)  - Combined LAI: 60 days observed + 10 days predicted
    !   current_doy    - Current day of year (1-365)
    !   current_state  - Current phenological state (STATE_DORMANT or STATE_GROWING)
    !
    ! Output:
    !   trigger_sos    - .true. if Start of Season should be triggered
    !   trigger_eos    - .true. if End of Season should be triggered
    !
    !-----------------------------------------------------------------------
      
      implicit none
      
      ! Arguments
      real(r8), intent(in)  :: lai_70day(N_TOTAL)
      integer,  intent(in)  :: current_doy
      integer,  intent(in)  :: current_state
      logical,  intent(out) :: trigger_sos
      logical,  intent(out) :: trigger_eos
      
      ! Local variables
      real(r8) :: deriv(N_TOTAL)
      integer  :: search_start, search_end
      integer  :: extremum_idx
      real(r8) :: extremum_val
      logical  :: found_extremum
      logical  :: in_sos_season, in_eos_season
      
      ! Initialize outputs
      trigger_sos = .false.
      trigger_eos = .false.
      
      ! Define search window: ±TOLERANCE days around TODAY_IDX
      search_start = TODAY_IDX - TOLERANCE   ! Day 55
      search_end   = TODAY_IDX + TOLERANCE   ! Day 65
      
      ! Compute first derivative of the 70-day LAI series
      call compute_first_derivative(lai_70day, N_TOTAL, deriv)
      
      ! Determine seasonal windows
      ! SOS season: from winter solstice to summer solstice (roughly DOY 1-171 in Northern Hemisphere)
      ! EOS season: from summer solstice to winter solstice (roughly DOY 172-365)
      in_sos_season = (current_doy >= 1) .and. (current_doy <= SUMMER_SOLSTICE)
      in_eos_season = (current_doy > SUMMER_SOLSTICE) .and. (current_doy <= 365)
      
      !-------------------------------------------------------------------
      ! Check for SOS (Start of Season)
      !-------------------------------------------------------------------
      if (in_sos_season .and. current_state == STATE_DORMANT) then
        
        ! Look for local MAXIMUM in derivative (peak growth rate = inflection point)
        call find_local_extremum(deriv, N_TOTAL, search_start, search_end, &
                                 .true., extremum_idx, extremum_val, found_extremum)
        
        ! Check if found and exceeds threshold
        ! Additional constraints:
        ! 1. Not too early in season (avoid Jan 1 false triggers)
        ! 2. LAI should be relatively low (< 2.0) at SOS
        ! 3. LAI should be increasing (check that we're coming from lower LAI)
        if (found_extremum .and. extremum_val > SOS_DERIV_THRESH .and. &
            current_doy > 30 .and. &
            lai_70day(TODAY_IDX) < 2.0_r8 .and. &
            lai_70day(TODAY_IDX) > lai_70day(TODAY_IDX-10)) then
          trigger_sos = .true.
        end if
        
      end if
      
      !-------------------------------------------------------------------
      ! Check for EOS (End of Season)
      !-------------------------------------------------------------------
      if (in_eos_season .and. current_state == STATE_GROWING) then
        
        ! Look for local MINIMUM in derivative (peak decline rate = inflection point)
        call find_local_extremum(deriv, N_TOTAL, search_start, search_end, &
                                 .false., extremum_idx, extremum_val, found_extremum)
        
        ! Check if found and exceeds threshold (note: threshold is negative)
        ! Additional constraints:
        ! 1. LAI should be reasonably high (> 2.5) at EOS
        ! 2. LAI should be decreasing (check that we're going toward lower LAI)
        if (found_extremum .and. extremum_val < EOS_DERIV_THRESH .and. &
            lai_70day(TODAY_IDX) > 2.5_r8 .and. &
            lai_70day(TODAY_IDX) > lai_70day(TODAY_IDX+5)) then
          trigger_eos = .true.
        end if
        
      end if
      
    end subroutine detect_phenology_transition

      !-----------------------------------------------------------------------
    subroutine get_sos(lai, n, sos_flag)
      real(8), dimension(n), intent(in) :: lai
      integer, intent(in) :: n
      integer, intent(out) :: sos_flag
      real(8), dimension(n) :: t, kv
      real(8) :: a, b, c, d, max_annual_lai, annual_lai_amp
      real(8) :: z, z1, z2, z3, p1n, p1d, p2n, p2d
      integer :: i
      
      ! Initialize variables
      max_annual_lai = 1
      annual_lai_amp = 0.5
      do i = 1, n
        t(i) = real(i - 1, 8)
      end do
      
      ! Initial parameter guesses
      a = -5.0
      b = 0.1
      c = maxval(lai)
      d = minval(lai)
      
      ! Calculate curvature rate of change
      do i = 1, size(t)
        z = exp(a + b * t(i))
        z1 = 1.0 - z
        z2 = 1.0 + z
        z3 = (b * c * z)**2
        p1n = 3.0 * z * z1 * z2**3 * (2.0 * z2**3 + b**2 * c**2 * z)
        p1d = (z2**4 + z3)**(2.5)
        p2n = z2**2 * (1.0 + 2.0 * z - 5.0 * z**2)
        p2d = (z2**4 + z3)**(1.5)
        kv(i) = b**3 * c * z * ((p1n / p1d) - (p2n / p2d))
      end do
      
      if (maxval(lai) > max_annual_lai * 0.3 .and. maxval(lai) - minval(lai) > 0.3 * annual_lai_amp) then
        sos_flag = 0
        do i = 2, n - 1
          if (kv(i) > kv(i - 1) .and. kv(i) > kv(i + 1)) then
            sos_flag = 1
            return
          end if
        end do
      else
        sos_flag = 0
      end if
      
    end subroutine get_sos

  !-----------------------------------------------------------------------
    function SeasonalDecidOnset( onset_gdd, onset_gddflag, soilt, dayofyear ) &
                       result( do_onset )

      ! !DESCRIPTION:
      ! Function to determine if seasonal deciduous leaf onset should happen.
      !
      ! !ARGUMENTS:
      real(r8), intent(INOUT) :: onset_gdd      ! onset growing degree days 
      real(r8), intent(INOUT) :: onset_gddflag  ! Onset freeze flag
      real(r8), intent(IN)    :: soilt          ! Soil temperature at specific level for this evaluation
      real(r8), intent(IN)    :: dayofyear      ! day of year
      logical :: do_onset                       ! Flag if onset should happen (return value)

      ! !LOCAL VARIABLES:
      real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
      real(r8):: crit_onset_gdd !critical onset growing degree-day sum
      real(r8):: crit_dayl      ! parameter
      real(r8):: annavg_t2m_patch        

      !-----------------------------------------------------------------------
      ! set constants
      annavg_t2m_patch = 15 ! annual average patch temperature (C)
      crit_dayl = 39300 ! seconds
      
      ! onset gdd sum from Biome-BGC, v4.1.2
      crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m_patch))
  
      ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
      if (dayofyear <= 171) then
        ws_flag = 1._r8
      else
        ws_flag = 0._r8
      end if
  
      do_onset = .false.
      ! Test to turn on growing degree-day sum, if off.
      ! switch on the growing degree day sum on the winter solstice
  
      if (onset_gddflag == 0._r8 .and. ws_flag == 1._r8) then
        onset_gddflag = 1._r8
        onset_gdd = 0._r8
      end if
  
      ! Test to turn off growing degree-day sum, if on.
      ! This test resets the growing degree day sum if it gets past
      ! the summer solstice without reaching the threshold value.
      ! In that case, it will take until the next winter solstice
      ! before the growing degree-day summation starts again.
  
      if (onset_gddflag == 1._r8 .and. ws_flag == 0._r8) then
        onset_gddflag = 0._r8
        onset_gdd = 0._r8
      end if
  
      ! if the gdd flag is set, and if the soil is above freezing
      ! then accumulate growing degree days for onset trigger
  
      if (onset_gddflag == 1.0_r8 .and. soilt > 273.15_r8) then
        onset_gdd = onset_gdd + (soilt-273.15_r8)
      end if

      ! set do_onset if critical growing degree-day sum is exceeded
      if (onset_gdd > crit_onset_gdd) then
        do_onset = .true.
      end if
  
    end function SeasonalDecidOnset

    !-----------------------------------------------------------------------
    subroutine next_ten_days(dayofyear, doy_arr)
      real(r8), intent(in)  :: dayofyear
      real(r8), intent(out) :: doy_arr(10)
      integer :: i
      real(r8) :: mod_doy
      mod_doy = real(365, r8)

      do i = 1, 10
        doy_arr(i) = mod(dayofyear + real(i, r8) - real(1, r8), mod_doy) + real(1, r8)
      end do
    end subroutine next_ten_days

    subroutine run_tft_model (the_tft_torch_model, ta_min, pr, sw, lai, doy_arr, out_data_tft)

      use   iso_c_binding,     only : c_float, c_int
      use   ftorch,            only : torch_model, torch_model_load, torch_model_forward, &
                                      torch_tensor, torch_tensor_from_array, torch_kCPU,  torch_delete
      implicit none
  
      ! Arguments
      character(len=*), intent(in) :: the_tft_torch_model
      real(r8),         intent(in) :: ta_min(:), pr(:), sw(:), lai(:), doy_arr(10)
      real(r8),         intent(out) :: out_data_tft(1,10)

      ! Local
      type(torch_model)                             :: model_pytorch
      type(torch_tensor), allocatable               :: inputs(:)
      type(torch_tensor), allocatable               :: outputs(:)
      real(c_float),     dimension(0)               :: empty_cat_raw
      real(c_float),      dimension(1,2)            :: static_num                           ! latitude, longitude
      integer(c_int),     dimension(1,0)            :: static_cat                           ! no static categorical features
      real(c_float),      dimension(1,60,8)         :: hist_num                             ! tmin, tmax, precip, rad, photoperiod, swvl1, doy, lai
      integer(c_int),     dimension(1,60,0)         :: hist_cat                             ! no historical categorical features
      real(c_float),      dimension(1,10,1)         :: fut_num                              ! doy
      integer(c_int),     dimension(1,10,0)         :: fut_cat                              ! no future categorical features
      real(c_float),      dimension(1,10,3)         :: out_quantiles                        ! output quantiles (0.1, 0.5, 0.9)
      integer(c_int),     allocatable               :: L2(:), L3(:)
      integer                                       :: n_in, n_out, i
      real(c_float), dimension(8) :: hmin, hmax
      real(c_float)              :: fmin, fmax
      integer                    :: j

      n_in  = 6                                                                             ! 6 input tensors
      n_out = 1                                                                             ! 1 output tensor (the 3‐quantile forecast)

      !---  Populate input data (first n_input days)
      static_num= reshape([ 1_c_float, 1_c_float ], [1,2])    ! TD: substitute with true lat/lon

      hist_num(1,:,1)= real(ta_min(1:60), c_float)                        ! TD: replace with tmin
      hist_num(1,:,2)= real(ta_min(1:60), c_float)                        ! TD: replace with tmax 
      hist_num(1,:,3)= real(pr(1:60), c_float)
      hist_num(1,:,4)= real(sw(1:60), c_float)
      ! TODO: Replace the following placeholder assignments with the correct variables for each feature
      hist_num(1,:,5)= real(ta_min(1:60), c_float)                        ! TODO: replace with photoperiod variable
      hist_num(1,:,6)= real(ta_min(1:60), c_float)                        ! TODO: replace with soil moisture variable
      hist_num(1,:,7)= real(ta_min(1:60), c_float)                        ! TODO: replace with day of year variable
      hist_num(1,:,8)= real(ta_min(1:60), c_float)                        ! TODO: replace with lai variable

      fut_num(1,:,1) = real(doy_arr(1:10), c_float)                   ! future 10 days of year

      static_cat    = reshape(empty_cat_raw, [1,0])                   ! no static categorical features
      hist_cat      = reshape(empty_cat_raw, [1,60,0])                ! no historical categorical features
      fut_cat       = reshape(empty_cat_raw, [1,10,0])                ! no future categorical features

      ! ---- normalize historical channels ----
      do j = 1, 8
        hmin(j) = minval( hist_num(1, :, j) )
        hmax(j) = maxval( hist_num(1, :, j) )
        if (hmax(j) > hmin(j)) then
          hist_num(1, :, j) = ( hist_num(1, :, j) - hmin(j) ) / (hmax(j) - hmin(j))
        else
          hist_num(1, :, j) = 0.0_c_float  ! or leave at 0 if flat
        end if
      end do

      ! ---- normalize future single channel ----
      fmin = minval( fut_num(1, :, 1) )
      fmax = maxval( fut_num(1, :, 1) )
      if (fmax > fmin) then
        fut_num(1, :, 1) = ( fut_num(1, :, 1) - fmin ) / (fmax - fmin)
      else
        fut_num(1, :, 1) = 0.0_c_float
      end if

      !===============
      ! load pytorch model
      call torch_model_load(model_pytorch, trim(the_tft_torch_model), torch_kCPU)

      !===============
      ! Allocate arrays of tensor handles
      allocate(inputs(n_in))
      allocate(outputs(n_out))

      allocate(L2(2));  L2  = [1_c_int,2_c_int]
      allocate(L3(3)); L3 = [1_c_int,2_c_int,3_c_int]
              
      !===============
      ! Wrap each Fortran array in a torch_tensor
      call torch_tensor_from_array(inputs(1), static_num, L2, torch_kCPU)
      call torch_tensor_from_array(inputs(2), static_cat, L2, torch_kCPU)
      call torch_tensor_from_array(inputs(3), hist_num,   L3, torch_kCPU)
      call torch_tensor_from_array(inputs(4), hist_cat,   L3, torch_kCPU)
      call torch_tensor_from_array(inputs(5), fut_num,    L3, torch_kCPU)
      call torch_tensor_from_array(inputs(6), fut_cat,    L3, torch_kCPU)

      !===============
      ! Wrap the output buffer
      call torch_tensor_from_array(outputs(1), out_quantiles, L3, torch_kCPU)
  
      !===============
      ! run pytorch model
      call torch_model_forward(model_pytorch, inputs, outputs)
      !===============
      ! extract median quantile (0.5)
      out_data_tft(1,:) = real(out_quantiles(1,:,2), kind=r8)
      !===============
      ! free pytorch model and tensors
      ! call torch_model_free(model_pytorch)                                       ! TO DO: confirm torch_model_free in ftorch
      do i = 1, n_in
        call torch_delete(inputs(i))
      end do
      call torch_delete(outputs(1))
      deallocate(inputs, outputs, L2, L3)

    end subroutine run_tft_model

    !-----------------------------------------------------------------------

    subroutine run_phenology_loop(lai_full, n_days, doy_full, &
                                 sos_doy, eos_doy, n_years)
      !
      ! Main driver that loops through the forcing data and detects 
      ! phenological transitions using the ML-based method.
      !
      ! This is a simplified version for testing that:
      !   - Uses observed LAI for both historical and "predicted" windows
      !   - In the full implementation, the 10-day predicted window would 
      !     come from the TFT model
      !
      ! Input:
      !   lai_full(n_days)  - Full LAI time series from forcing data
      !   n_days            - Total number of days in forcing data
      !   doy_full(n_days)  - Day of year for each timestep
      !
      ! Output:
      !   sos_doy(n_years)  - Detected SOS day of year for each year
      !   eos_doy(n_years)  - Detected EOS day of year for each year
      !   n_years           - Number of complete years processed
      !
      !-----------------------------------------------------------------------
      
      implicit none
      
      ! Arguments
      integer,  intent(in)  :: n_days
      real(r8), intent(in)  :: lai_full(n_days)
      real(r8), intent(in)  :: doy_full(n_days)
      integer,  intent(out) :: n_years
      integer,  intent(out) :: sos_doy(:)
      integer,  intent(out) :: eos_doy(:)
      
      ! Local variables
      real(r8) :: lai_70day(N_TOTAL)
      integer  :: current_state
      integer  :: current_doy
      integer  :: day_idx
      integer  :: year_idx
      integer  :: last_doy
      logical  :: trigger_sos, trigger_eos
      
      ! Initialize
      current_state = STATE_DORMANT   ! Start in dormant state
      year_idx = 1
      n_years = 0
      sos_doy = 0
      eos_doy = 0
      last_doy = 0
      
      ! Loop through days, starting at day 61 (need 60 days of history)
      do day_idx = N_HIST + 1, n_days - N_PRED
        
        current_doy = int(doy_full(day_idx))
        
        ! Detect year transition (DOY goes from 365/366 back to 1)
        if (current_doy < last_doy) then
          year_idx = year_idx + 1
          ! Reset state at start of new year if still growing (shouldn't happen normally)
          if (current_state == STATE_GROWING) then
            current_state = STATE_DORMANT
          end if
        end if
        last_doy = current_doy
        
        ! Check array bounds for output
        if (year_idx > size(sos_doy)) exit
        
        ! Extract 70-day window
        ! Historical: days (day_idx - 60) to (day_idx - 1)
        ! "Today" + predicted: days day_idx to (day_idx + 9)
        lai_70day(1:N_HIST) = lai_full(day_idx - N_HIST : day_idx - 1)
        lai_70day(N_HIST+1:N_TOTAL) = lai_full(day_idx : day_idx + N_PRED - 1)
        
        ! Detect transitions
        call detect_phenology_transition(lai_70day, current_doy, current_state, &
                                          trigger_sos, trigger_eos)
        
        ! Handle SOS trigger
        if (trigger_sos) then
          sos_doy(year_idx) = current_doy
          current_state = STATE_GROWING
          print *, "  SOS triggered on DOY ", current_doy, " (year ", year_idx, ")"
        end if
        
        ! Handle EOS trigger
        if (trigger_eos) then
          eos_doy(year_idx) = current_doy
          current_state = STATE_DORMANT
          print *, "  EOS triggered on DOY ", current_doy, " (year ", year_idx, ")"
          n_years = year_idx  ! Mark this year as complete
        end if
        
      end do
      
      ! If we ended mid-year, count that year if SOS was detected
      if (sos_doy(year_idx) > 0 .and. n_years < year_idx) then
        n_years = year_idx
      end if
      
    end subroutine run_phenology_loop
  
end program MLPhenology
