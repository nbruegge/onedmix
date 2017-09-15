module onedmix_vmix_mytke
  use onedmix_variables
  use onedmix_eos
  use onedmix_io
  implicit none

  ! namelist parameters
  real*8                :: alpha_tke, c_eps, cu, cd
  real*8                :: Lmix_min, kappaM_min, kappaM_max, tke_min

  ! TKE diagnostics
  real*8, dimension(:), allocatable :: &
                                       Etke,          &
                                       Gimp_tke,      &
                                       Gexp_tke,      &
                                       TEtke_dif,     &
                                       TEtke_dis,     &
                                       TEtke_spr,     &
                                       TEtke_bpr,     &
                                       TEtke_tau,     &
                                       TEtke_tot
  
  contains
!-------------------------------------------------------------------------------- 
  subroutine setup_vmix_mytke
    namelist/tke/            alpha_tke, c_eps, cu, cd, &
                             Lmix_min, kappaM_min, kappaM_max, tke_min
    ! read namelist or take standard parameters
    if (.true.) then
      open(fid, file="./onedmix.nl", status="old", action='read')
      read(fid, nml=tke)
      close(fid)
    end if
    !write(*,*) 'alpha_tke = ', alpha_tke
    !write(*,*) 'c_eps = ', c_eps
    !write(*,*) 'cu = ', cu
    !write(*,*) 'cd = ', cd
    !write(*,*) 'Lmix_min = ', Lmix_min
    !write(*,*) 'kappaM_min = ', kappaM_min
    !write(*,*) 'kappaM_max = ', kappaM_max
    !write(*,*) 'tke_min = ', tke_min

    allocate( Etke(1:nz+1) ); Etke=0.0
    allocate( Gimp_tke(1:nz+1) ); Gimp_tke=0.0
    allocate( Gexp_tke(1:nz+1) ); Gexp_tke=0.0

    allocate( TEtke_dif(1:nz+1) ); TEtke_dif=0.0
    allocate( TEtke_dis(1:nz+1) ); TEtke_dis=0.0
    allocate( TEtke_spr(1:nz+1) ); TEtke_spr=0.0
    allocate( TEtke_bpr(1:nz+1) ); TEtke_bpr=0.0
    allocate( TEtke_tau(1:nz+1) ); TEtke_tau=0.0
    allocate( TEtke_tot(1:nz+1) ); TEtke_tot=0.0

  end subroutine setup_vmix_mytke

!-------------------------------------------------------------------------------- 
  subroutine calc_vmix_mytke
    real*8, dimension(nz+1) :: Gimp_tke_new, Gexp_tke_new, Etke_exp
    real*8, dimension(nz+1) :: adiff, bdiff, cdiff, ldiag, mdiag, udiag
    real*8, dimension(nz+1) :: ke, Pr, cb
    real*8, dimension(nz+1) :: delta, Lmix, sqrtEtke
    real*8, dimension(nz+1) :: K_diss_v, P_diss_v
    real*8, dimension(nz+1) :: old_Etke
    real*8                :: forc_rho_surf, forc_tke_surf
    integer :: k

    adiff = 0.0
    bdiff = 0.0
    cdiff = 0.0
    ldiag = 0.0
    mdiag = 0.0
    udiag = 0.0
    Gimp_tke_new = 0.0
    Gexp_tke_new = 0.0
    Etke_exp = 0.0

    ! FIXME: Derive this!
    forc_rho_surf = 0.0
    !forc_tke_surf = 2e-5
    !taux = 2.4683998517220244E-004
    !tauy = 1.2652587216272979E-004
    !taux_act = 1.0550815445478266E-004
    !tauy_act = 9.9133200202586159E-006
    old_Etke = Etke

    sqrtEtke = sqrt(max(0d0,Etke))
    Lmix = sqrt(2.0)*sqrtEtke/sqrt(max(1d-12, N2))
    ! !!! Do mixing length restriction stuff
    ! constrain mixing length scale as in MITgcm/OPA code (from pyOM)
    do k=2,nz+1
      Lmix(k) = MIN(Lmix(k), Lmix(k-1)+dzt(k-1) )
    end do 
    Lmix(1) = MIN( Lmix(1), Lmix_min+dzt(1))
    do k=nz,1,-1
      Lmix(k) = MIN(Lmix(k), Lmix(k+1)+dzt(k))
    end do
    Lmix = max(Lmix, Lmix_min)

    ! Pr=6.6*Ri with 1<Pr<10
    Pr = max(1.0, min(10.0,6.6*Ri) )
    Av = 0.
    kv = 0.
    !Av = Av + cu * sqrt(Etke) * Lmix
    !kv = kv + cb * sqrt(Etke) * Lmix
    Av = min(kappaM_max, cu*sqrtEtke*Lmix )
    kv = Av / Pr
    !Av = Av + Avb
    !kv = kv + kvb

    ! diffusion of Etke
    ke = alpha_tke*Av
    do k=2,nz+1
      delta(k) = 0.5*(ke(k)+ke(k-1))/dzw(k-1)
    end do
    delta(1) = 0.0

    do k=2,nz+1
      !adiff(k) = ke(k)/(dzt(k)*dzw(k))
      adiff(k) = delta(k)/dzt(k)
    end do
    adiff(1) = 0.0 ! since ke(1)=0.
    do k=2,nz
      !bdiff(k) = ke(k)/(dzt(k)*dzw(k)) + ke(k+1)/(dzt(k)*dzw(k+1))
      bdiff(k) = delta(k)/dzt(k) + delta(k+1)/dzt(k)
    end do
    !bdiff(1)  = ke(2)/(dzw(2)*dzt(2))
    !bdiff(nz) = ke(nz)/(dzw(nz)*dzt(nz))
    bdiff(1)  = delta(2)/dzt(1)
    bdiff(nz+1) = delta(nz+1)/dzt(nz+1)
    do k=1,nz
      !cdiff(k) = ke(k+1)/(dzt(k)*dzw(k+1))
      cdiff(k) = delta(k+1)/dzt(k)
    end do
    cdiff(nz+1) = 0.0 

    ! implicit part
    do k=2,nz
      Gimp_tke_new(k) = adiff(k)*Etke(k-1) - bdiff(k)*Etke(k) + cdiff(k)*Etke(k+1)
    end do
    Gimp_tke_new(1)    = cdiff(1)*Etke(2)     - bdiff(1)*Etke(1)
    Gimp_tke_new(nz+1) = adiff(nz+1)*Etke(nz) - bdiff(nz+1)*Etke(nz+1)
    TEtke_dif = Gimp_tke_new
    ! add Etke dissipation
    Gimp_tke_new = Gimp_tke_new - c_eps*sqrtEtke**3/Lmix
    TEtke_dis = -c_eps*sqrtEtke**3/Lmix

    ! explicit part
    ! add forcing
    P_diss_v   = kv*N2
    K_diss_v   = Av*S2
    P_diss_v(1) = -forc_rho_surf*grav/rho0
    TEtke_spr = K_diss_v
    TEtke_bpr = -P_diss_v
    Gexp_tke_new = K_diss_v-P_diss_v
    forc_tke_surf = cd * sqrt( (taux_act**2 + tauy_act**2)**(3./2.)  )
    Gexp_tke_new(1) = Gexp_tke_new(1) + forc_tke_surf/dzt(1)
    TEtke_tau(1) = forc_tke_surf/dzt(1)

    ! Adams-Bashforth time stepping
    !write(*,*) 'Etke_exp = ', Etke_exp
    !write(*,*) 'Etke = ', Etke
    !write(*,*) 'Gimp_tke_new = ', Gimp_tke_new
    !write(*,*) 'Gimp_tke = ', Gimp_tke
    !write(*,*) 'Gexp_tke_new = ', Gexp_tke_new
    !write(*,*) 'Gexp_tke = ', Gexp_tke
    Etke_exp = Etke+dt*(1.d0-dimp)*((1.5+epsab)*Gimp_tke_new-(0.5+epsab)*Gimp_tke) &
                   +dt*            ((1.5+epsab)*Gexp_tke_new-(0.5+epsab)*Gexp_tke)
    Gimp_tke = Gimp_tke_new
    Gexp_tke = Gexp_tke_new

    ! implicit part
    ! FIXME: Add dissipation and forcing to mdiag
    udiag = -dt*dimp*cdiff
    mdiag = (1+dt*dimp*(bdiff + c_eps*sqrtEtke/Lmix))
    ldiag = -dt*dimp*adiff
    call solve_tridiag(ldiag, mdiag, udiag, Etke_exp, Etke, nz+1)
 
    ! set bounding values of Etke
    if (Etke(1) < 0.0) then
      Etke(1) = 0.0
    end if
    Etke(2:nz+1) = max(Etke(2:nz+1), tke_min)
    TEtke_tot = (Etke-old_Etke)/dt

    tstep_count = tstep_count+1
    if (.false.) then
      write(*,*) "================================================================================"  
      write(*,*) 'tstep_count = ', tstep_count
      !write(*,*) 'uvel = ', uvel
      !write(*,*) 'vvel = ', vvel
      !write(*,*) 'S2 = ', S2
      !stop
      !write(*,*) 'temp = ', temp
      !write(*,*) 'salt = ', salt
      write(*,*) 'N2 = ', N2
      write(*,*) 'S2 = ', S2
      !write(*,*) 'dzw = ', dzw
      !write(*,*) 'dzt = ', dzt
      !write(*,*) 'Pr  = ', Pr
      write(*,*) 'Lmix = ', Lmix
      write(*,*) 'sqrtEtke = ', sqrtEtke
      write(*,*) 'Av = ', Av
      write(*,*) 'delta*dt = ', delta*dt
      write(*,*) 'Gexp_tke_new = ', Gexp_tke_new
      !write(*,*) 'forc_tke_surf = ', forc_tke_surf
      write(*,*) 'ldiag = ', ldiag
      write(*,*) 'mdiag = ', mdiag
      write(*,*) 'udiag = ', udiag
      write(*,*) 'Etke_exp = ', Etke_exp
      write(*,*) 'Etke = ', Etke
      !write(*,*) 'kv = ', kv
      write(*,*) "================================================================================"  
    if ( tstep_count==1 ) then
      stop
    end if
    end if

  end subroutine calc_vmix_mytke

!-------------------------------------------------------------------------------- 
  subroutine write_snap_mytke
    character(len=20)     :: fprfx

    fprfx = 'onedmix_state       '
    call save_variable(fprfx, Etke, 'Etke', 'turbulent kinetic energy', &
                       iostep, nz+1, 'm^2 / s^2')

    call save_variable(fprfx, TEtke_dif, 'TEtke_dif', &
                       'tke tendency by diffusion', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, TEtke_dis, 'TEtke_dis', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, TEtke_spr, 'TEtke_spr', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, TEtke_bpr, 'TEtke_bpr', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, TEtke_tau, 'TEtke_tau', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, TEtke_tot, 'TEtke_tot', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')
  end subroutine write_snap_mytke

end module onedmix_vmix_mytke
