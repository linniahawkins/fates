program MLPhenology

  use               FatesConstantsMod, only : r8 => fates_r8
  use              FatesArgumentUtils, only : command_line_arg
  use                FatesTestFireMod, only : ReadDatmData
  !use               ftorch,            only : torch_tensor, torch_tensor_from_array, torch_kCPU
  !use               ftorch,            only : torch_model, torch_model_load, torch_model_forward, &
  !                                            torch_tensor, torch_tensor_from_array, torch_kCPU,  torch_delete  

  implicit none

  !'/glade/u/home/linnia/MLphenology/US-MMS_forcing_2001.nc'

  ! define ML model
  !character(len=256) :: nn_torch_model = "/glade/u/home/linnia/FTorch_example/saved_simplenet_model_v2.pt"
  !character(len=256) :: lstm_torch_model = "/glade/u/home/ayal/phenology-ml-clm/models/example_LSTM_model_v1.pt"
  !type(torch_model) :: model_pytorch                                    

  ! setup Fortran data structures  
  !type(torch_tensor), dimension(1)         :: in_tensor, out_tensor
  !integer                                  :: in_layout(1) = [1]
  !integer                                  :: out_layout(1) = [1] 

  character(len=:),                  allocatable :: datm_file            ! input DATM 
  real(r8),                          allocatable :: temp_degC(:)         ! daily air temperature [degC]
  real(r8),                          allocatable :: precip(:)            ! daily precipitation [mm]
  real(r8),                          allocatable :: sw(:)                ! daily shortwave radiation (W/m2)
  real(r8),                          allocatable :: lai(:)              ! daily LAI (m2/m2)

  real(r8),              dimension(60,4), target :: in_data
  real(r8),                 dimension(5), target :: out_data

  ! allocate arrays
  allocate(temp_degC(n_days))
  allocate(precip(n_days))
  allocate(sw(n_days))
  allocate(lai(n_days))

  ! read in DATM data
  datm_file = command_line_arg(1)
  call ReadDatmData(datm_file, temp_degC, precip, rh, wind)
 
  !===============
  ! load model
  !call torch_model_load(model_pytorch, trim(nn_torch_model), torch_kCPU)
  
  ! set input data 
  in_data = 0.5_r8
  call random_number(in_data(:,1))   ! fill first column with random numbers in [0,1)

  print *, in_data(1:5, :)
  
  !print *, in_data
  !print *, 'Shape of input_data:', size(in_data, 1), 'rows ×', size(in_data, 2), 'columns'
  

  !call torch_tensor_from_array(in_tensor(1), in_data, in_layout, torch_kCPU)
  !in_tensor(1) = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
  !print *, "Tensor shape: ", shape(in_tensor)
  
  
  !call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
  !call torch_model_forward(model_pytorch, in_tensor, out_tensor)

  !call torch_delete(in_tensors(1))
  !call torch_delete(out_tensors(1)) 

  print *, "Hello Phenology"
    
end program MLPhenology
