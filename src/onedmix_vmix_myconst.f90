module onedmix_vmix_myconst
  use onedmix_variables
  use onedmix_eos
  implicit none
  
  ! namelist parameters
  logical :: convective_adjustment
  real*8 :: kv_cadjust

  contains

!-------------------------------------------------------------------------------- 
  subroutine setup_vmix_myconst
    convective_adjustment = .false.
    kv_cadjust = 0.1
  end subroutine setup_vmix_myconst

!-------------------------------------------------------------------------------- 
  subroutine calc_vmix_myconst
    integer :: k
    Av = Avb
    kv = kvb

    ! --- convective adjustment
    if ( convective_adjustment ) then
      do k=1,nz
        if (N2(k)<0.0) then
          kv = kv_cadjust
        end if
      end do
    end if
  end subroutine calc_vmix_myconst

!-------------------------------------------------------------------------------- 
  subroutine write_snap_myconst
  end subroutine write_snap_myconst

end module onedmix_vmix_myconst
