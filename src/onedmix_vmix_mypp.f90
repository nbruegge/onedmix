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
    real*8 :: pint, dens_km1, dens_k
    integer :: k

    ! --- initialize/reset values
    N2=0.0
    S2=0.0
    uz=0.0
    vz=0.0
    Ri=0.0

    do k=2,nz
      !pint = 0.5*(pbcl(k-1)+pbcl(k))
      !pint = grav*rho0*zu(k)
      ! as in MPIOM beleg.f90 (l. 79)
      pint = 0.0001 * rho0 * zu(k)
      call potrho(temp(k-1), salt(k-1), pint, dens_km1) 
      call potrho(temp(k),   salt(k),   pint, dens_k) 
      N2(k) = -grav/rho0*(dens_km1-dens_k)/dzt(k)
      uz(k) = (uvel(k-1)-uvel(k))/dzt(k)
      vz(k) = (vvel(k-1)-vvel(k))/dzt(k)
    end do
    N2(1) = 0.0
    uz(1) = 0.0
    vz(1) = 0.0
    N2(nz+1) = 0.0
    uz(nz+1) = 0.0
    vz(nz+1) = 0.0
    S2 = uz**2+vz**2

    Ri = N2/(S2+1e-30)

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
