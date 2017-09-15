module onedmix_vmix_mypp
  use onedmix_variables
  use onedmix_eos
  implicit none
  
  ! namelist parameters
  real*8 :: Av0, kv0, alpha, nAv
  logical :: convective_adjustment
  real*8 :: kv_cadjust

  contains

!-------------------------------------------------------------------------------- 
  subroutine setup_vmix_mypp
    ! set default values
    Av0 = 5e-3
    kv0 = 5e-3
    alpha = 5.0
    nAv = 2.0

    convective_adjustment = .false.
    kv_cadjust = 10.0
  end subroutine setup_vmix_mypp
  
!-------------------------------------------------------------------------------- 
  subroutine calc_vmix_mypp
    integer :: k

    ! --- initialize/reset values
    Av=0.0
    kv=0.0
    Av = Av + Av0 / (1+alpha*Ri)**nAv
    kv = kv + kv0 / (1+alpha*Ri)
    Av = Av + Avb
    kv = kv + kvb

    ! --- convective adjustment
    if ( convective_adjustment ) then
      do k=1,nz
        if (N2(k)<0.0) then
          kv = 1e2
        end if
      end do
    end if
  end subroutine calc_vmix_mypp

!-------------------------------------------------------------------------------- 
  subroutine write_snap_mypp
  end subroutine write_snap_mypp

end module onedmix_vmix_mypp
