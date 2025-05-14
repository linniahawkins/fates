program MLPhenology

  use               FatesConstantsMod, only : r8 => fates_r8
  use              FatesArgumentUtils, only : command_line_arg
  use             FatesUnitTestIOMod,  only : OpenNCFile, GetVar, CloseNCFile, RegisterNCDims

  implicit none

  ! define ML phenoogy pytorch model
  !character(len=256) :: the_torch_model = "/glade/u/home/linnia/FTorch_example/constant_model.pt"
  character(len=256) :: the_torch_model = "/glade/u/home/ayal/phenology-ml-clm/models/example_LSTM_model_v1.pt"
  
  real(8), dimension(10) :: dummy_lai
  integer :: sos_flag, n

  character(len=:),                  allocatable :: datm_file            ! input DATM 
  real(r8),                          allocatable :: ta(:)         ! daily air temperature [degC]
  real(r8),                          allocatable :: pr(:)            ! daily precipitation [mm]
  real(r8),                          allocatable :: sw(:)                ! daily shortwave radiation (W/m2)
  real(r8),                          allocatable :: lai(:)              ! daily LAI (m2/m2)

  real(r8)                                       :: out_data(5)       ! output from the lstm model (lai)
  
  real(r8)                                       :: soilt            ! soil temperature at 12cm
  real(r8)                                       :: doy ! day of year (used to identify solstace) 
  real(r8)                                       :: onset_gdd      ! onset growing degree days 
  real(r8)                                       :: onset_gddflag  ! Onset freeze flag
  logical                                        :: do_onset       ! Flag if onset should happen

  
  ! Load forcing data
  datm_file = command_line_arg(1) ! one year of daily ta, pr, sw, lai
  call load_met_forcing(datm_file, ta, pr, sw, lai)

  ! ======================================
  ! test CLM SeasonalDecidOnset function
  doy = 1.0_r8
  onset_gdd = 0.0_r8
  onset_gddflag = 1.0_r8

  soilt = ta(doy)-10.0_r8

  do_onset = SeasonalDecidOnset( onset_gdd, onset_gddflag, soilt, doy )
  print *, "onset_gdd: ", onset_gdd
  print *, "onset_gddflag: ", onset_gddflag

  ! ====================================
  ! test SOS
  n = 10
  dummy_lai = (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0/)
  call get_sos(dummy_lai, n, sos_flag)
  print *, "start of season Flag:", sos_flag

  ! ========================================
  ! test lstm
  call run_pytorch_model(the_torch_model, ta, pr, sw, lai, doy, out_data)
  print *, out_data

  contains
  
    !-----------------------------------------------------------------------
    subroutine load_met_forcing ( datm_file, ta, pr, sw, lai)
      ! 
      use FatesConstantsMod, only: r8 => fates_r8
      use FatesUnitTestIOMod, only: OpenNCFile, GetVar, CloseNCFile
    
      implicit none
    
      ! Arguments
      character(len=*), intent(in) :: datm_file
      real(r8), allocatable, intent(out) :: ta(:), pr(:), sw(:), lai(:)
    
      ! Local
      integer :: ncid
    
      ! Allocate arrays
      allocate(ta(365), pr(365), sw(365), lai(365))
    
      ! Open and read
      call OpenNCFile(trim(datm_file), ncid, 'read')
      
      call GetVar(ncid, 'ta', ta)
      call GetVar(ncid, 'pr', pr)
      call GetVar(ncid, 'sw', sw)
      call GetVar(ncid, 'lai', lai)
    
      call CloseNCFile(ncid)
    
    end subroutine load_met_forcing
    
    ! ----------------------------------------------------------------------

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
    function SeasonalDecidOnset( onset_gdd, onset_gddflag, soilt, doy ) &
                       result( do_onset )

        ! !DESCRIPTION:
        ! Function to determine if seasonal deciduous leaf onset should happen.
        !
        ! !ARGUMENTS:
        real(r8), intent(INOUT) :: onset_gdd      ! onset growing degree days 
        real(r8), intent(INOUT) :: onset_gddflag  ! Onset freeze flag
        real(r8), intent(IN)    :: soilt          ! Soil temperature at specific level for this evaluation
        real(r8), intent(IN)    :: doy            ! day of year
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
        if (doy <= 171) then
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
    subroutine run_pytorch_model (the_torch_model, ta, pr, sw, lai, doy, out_data)

        use   FatesConstantsMod, only : r8 => fates_r8
        use   ftorch,            only : torch_model, torch_model_load, torch_model_forward, &
                                        torch_tensor, torch_tensor_from_array, torch_kCPU,  torch_delete  
        
        implicit none
    
        ! Arguments
        character(len=*), intent(in) :: the_torch_model
        real(r8),         intent(in) :: ta(:), pr(:), sw(:), lai(:)
        real(r8),         intent(in) :: doy            ! day of year
        real(r8),        intent(out) :: out_data(5)
    
        ! Local
        type(torch_model) :: model_pytorch
        type(torch_tensor), dimension(1)         :: in_tensor, out_tensor
        integer                                  :: in_layout(2) = [60, 4]
        integer                                  :: out_layout(1) = [5]
        real(r8),        dimension(60,4), target :: in_data

        ! Populate input data (first n_input days)
        in_data(:,1) = lai(1:60)
        in_data(:,2) = ta(1:60)
        in_data(:,3) = pr(1:60)
        in_data(:,4) = sw(1:60)
    
        !===============
        ! load pytorch model
        
        call torch_model_load(model_pytorch, trim(the_torch_model), torch_kCPU)
        
        !===============
        ! run pytorch model
        
        call torch_tensor_from_array(in_tensor(1), in_data, in_layout, torch_kCPU)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        call torch_model_forward(model_pytorch, in_tensor, out_tensor)
        
        call torch_delete(in_tensor(1))
        call torch_delete(out_tensor(1)) 
    
    end subroutine run_pytorch_model
        
end program MLPhenology
