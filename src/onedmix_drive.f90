program OneDmix 
  use OneDmix_setup
  use OneDmix_timeloop
  !use diffusion
  implicit none
  
  ! --- setup all model variables and forcing fields
  ! (OneDmix_setup/setup_OneDmix)
  call setup_OneDmix()

  ! --- main time loop
  ! (OneDmix_timeloop/timeloop)
  call timeloop()

  ! --- model run is completed
  write(*,*) 'All done!'
  write(*,*) "================================================================================"

end program OneDmix

