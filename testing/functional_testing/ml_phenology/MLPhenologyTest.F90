program MLPhenology

  use               FatesConstantsMod, only : r8 => fates_r8
  use               ftorch,            only : torch_tensor_from_array, torch_tensor, torch_kCPU
  
  implicit none 
    
  type(torch_tensor), dimension(1)         :: in_tensor
  integer                                  :: in_layout(1) = [1]
  real(r8),           dimension(1), target :: in_data
  real(8), dimension(10) :: lai
  integer :: sos_flag, n


  in_data = [1.0_r8]
  
  print *, "Hello Phenology"
    
  in_tensor(1) = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
  
  print *, "Tensor shape: ", shape(in_tensor)
  n = 10
  lai = (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0/)
  call get_sos(lai, n, sos_flag)
  print *, "SOS Flag:", sos_flag

  contains 
  ! Define the function to calculate the start of season (SOS)
  ! This function takes the LAI array and its size as input
  ! and returns the SOS flag (1 for SOS detected, 0 otherwise)
  ! The function uses a simple thresholding method to determine the SOS

  ! based on the maximum LAI value and its change over time
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
  
  
end program MLPhenology
