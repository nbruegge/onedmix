program onedmix 
  use onedmix_setup
  use onedmix_timeloop
  !use diffusion
  implicit none
  
  ! --- setup all model variables and forcing fields
  ! (onedmix_setup/setup_onedmix)
  call setup_onedmix()

  ! --- main time loop
  ! (onedmix_timeloop/timeloop)
  call timeloop()

  ! --- model run is completed
  write(*,*) 'All done!'
  write(*,*) "================================================================================"

end program onedmix

