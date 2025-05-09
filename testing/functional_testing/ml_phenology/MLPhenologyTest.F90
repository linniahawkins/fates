program MLPhenology

  use               FatesConstantsMod, only : r8 => fates_r8
  use              FatesArgumentUtils, only : command_line_arg
  use             FatesUnitTestIOMod,  only : OpenNCFile, GetVar, CloseNCFile, RegisterNCDims
  use                ftorch
  !use               ftorch,            only : torch_tensor, torch_tensor_from_array, torch_kCPU
  !use               ftorch,            only : torch_model, torch_model_load, torch_model_forward, &
  !                                            torch_tensor, torch_tensor_from_array, torch_kCPU,  torch_delete  

  implicit none


  ! define ML model
  character(len=256) :: cb_torch_model = "/glade/u/home/linnia/FTorch_example/constant_model.pt"
  !character(len=256) :: lstm_torch_model = "/glade/u/home/ayal/phenology-ml-clm/models/example_LSTM_model_v1.pt"
  type(torch_model) :: model_pytorch                                    

  ! setup Fortran data structures  
  type(torch_tensor), dimension(1)         :: in_tensor, out_tensor
  integer                                  :: in_layout(1) = [1]
  integer                                  :: out_layout(1) = [1] 
  integer,                       parameter :: n_days = 365
  integer                                  :: ncid ! netcdf file unit number

  character(len=:),                  allocatable :: datm_file            ! input DATM 
  real(r8),                          allocatable :: ta(:)         ! daily air temperature [degC]
  real(r8),                          allocatable :: pr(:)            ! daily precipitation [mm]
  real(r8),                          allocatable :: sw(:)                ! daily shortwave radiation (W/m2)
  real(r8),                          allocatable :: lai(:)              ! daily LAI (m2/m2)

  !real(r8)                                       :: in_data(1) = 1.0
  real(r8),              dimension(60,4), target :: in_data
  real(r8),                 dimension(1), target :: out_data

  ! allocate arrays
  allocate(ta(n_days))
  allocate(pr(n_days))
  allocate(sw(n_days))
  allocate(lai(n_days))

  !===============
  ! load meteorological forcing data

  ! open file
  datm_file = command_line_arg(1) ! one year of daily ta, pr, sw, lai
  call OpenNCFile(trim(datm_file), ncid, 'read')
      
  ! read in data
  call GetVar(ncid, 'ta', ta)
  call GetVar(ncid, 'pr', pr)
  call GetVar(ncid, 'sw', sw)
  call GetVar(ncid, 'lai', lai)   

  in_data(:,1)= lai(1:60)
  in_data(:,2)= ta(1:60)
  in_data(:,3)= pr(1:60)
  in_data(:,4)= sw(1:60)

  ! normalize input data for pytorch model

  print *, in_data(1,:)

  ! close file
  call CloseNCFile(ncid)

  !===============
  ! load pytorch model

  call torch_model_load(model_pytorch, trim(cb_torch_model), torch_kCPU)
  
  ! set input data 
  !in_data = 0.5_r8
  !call random_number(in_data(:,1))   ! fill first column with random numbers in [0,1)

  !print *, in_data(1:5, :)
  
  !print *, in_data
  !print *, 'Shape of input_data:', size(in_data, 1), 'rows ×', size(in_data, 2), 'columns'
  

  !===============
  ! run pytorch model

  !call torch_tensor_from_array(in_tensor(1), in_data, in_layout, torch_kCPU)
  !in_tensor(1) = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
  !print *, "Tensor shape: ", shape(in_tensor)
  !call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
  !call torch_model_forward(model_pytorch, in_tensor, out_tensor)

  !call torch_delete(in_tensor(1))
  !call torch_delete(out_tensor(1)) 

  print *, "Hello Phenology"
    
end program MLPhenology
