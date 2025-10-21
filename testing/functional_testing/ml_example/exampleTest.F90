program exampleTest

  use               FatesConstantsMod, only : r8 => fates_r8
  use              FatesArgumentUtils, only : command_line_arg
  use             FatesUnitTestIOMod,  only : OpenNCFile, GetVar, CloseNCFile, RegisterNCDims

  implicit none

  ! define ML phenology pytorch model
  character(len=256) :: the_torch_model = "/glade/u/home/linnia/MLphenology/models/example_LSTM_model_lh.pt"
  
  real(8), dimension(10) :: dummy_lai
  integer :: n

  character(len=:),                  allocatable :: datm_file            ! input DATM 
  real(r8),                          allocatable :: ta_min(:)         ! daily min air temperature [degC]
  real(r8),                          allocatable :: ta_max(:)         ! daily max air temperature [degC]
  real(r8),                          allocatable :: pr(:)            ! daily precipitation [mm]
  real(r8),                          allocatable :: sw(:)                ! daily shortwave radiation (W/m2)
  real(r8),                          allocatable :: lai(:)              ! daily LAI (m2/m2)
  real(r8),                          allocatable :: soilm(:)              ! daily soil moisture at layer 3 (kg/m2)
  real(r8),                          allocatable :: doy(:)              ! day of year
  real(r8),                          allocatable :: photo(:)              ! daily photoperiod (seconds)

  real(r8)                                       :: out_data(1,5)       ! output from the lstm model (lai)

  real(r8)                                       :: dayofyear ! day of year 
  real(r8)                                       :: soilt            ! soil temperature at 12cm 
  real(r8)                                       :: onset_gdd      ! onset growing degree days 
  real(r8)                                       :: onset_gddflag  ! Onset freeze flag
  logical                                        :: do_onset       ! Flag if onset should happen

  ! ======================================
  ! Load forcing data
  ! the path to this file is set in /src/fates/testing/functional_tests.cfg
  datm_file = command_line_arg(1) ! one year of daily ta, pr, sw, lai 
  call load_met_forcing(datm_file, ta_min, ta_max, pr, sw, lai, soilm, doy, photo)
  
  ! ========================================
  ! test lstm
  call run_pytorch_model(the_torch_model, ta_min, pr, sw, lai, dayofyear, out_data)
  print *, "predicted LAI:", out_data

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
    !-----------------------------------------------------------------------
    subroutine run_pytorch_model (the_torch_model, ta_min, pr, sw, lai, dayofyear, out_data)

        use   iso_c_binding,     only : c_float, c_int
        use   ftorch,            only : torch_model, torch_model_load, torch_model_forward, &
                                        torch_tensor, torch_tensor_from_array, torch_kCPU,  torch_delete  
        
        implicit none
    
        ! Arguments
        character(len=*), intent(in) :: the_torch_model
        real(r8),         intent(in) :: ta_min(:), pr(:), sw(:), lai(:)
        real(r8),         intent(in) :: dayofyear            ! day of year
        real(r8),        intent(out) :: out_data(1,5)
    
        ! Local
        type(torch_model) :: model_pytorch
        type(torch_tensor), dimension(1)         :: in_tensor, out_tensor
        integer(c_int)                                  :: in_layout(3) = [1,2,3]
        integer(c_int)                                  :: out_layout(2) = [1,2]
        real(c_float),        dimension(1,60,4), target :: in_data

        ! Populate input data (first n_input days)
        in_data(1,:,1) = real(lai(1:60), c_float)
        in_data(1,:,2) = real(ta_min(1:60), c_float)
        in_data(1,:,3) = real(pr(1:60), c_float)
        in_data(1,:,4) = real(sw(1:60), c_float)
    
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
        
end program exampleTest
