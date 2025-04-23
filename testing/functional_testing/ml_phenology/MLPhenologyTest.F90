program MLPhenology

  use               FatesConstantsMod, only : r8 => fates_r8
  use               ftorch,            only : torch_tensor_from_array, torch_tensor, torch_kCPU
  
  implicit none 
    
  type(torch_tensor), dimension(1)         :: in_tensor
  integer                                  :: in_layout(1) = [1]
  real(r8),           dimension(1), target :: in_data
  
  in_data = [1.0_r8]
  
  print *, "Hello Phenology"
    
  in_tensor(1) = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
  
  print *, "Tensor shape: ", shape(in_tensor)

end program MLPhenology