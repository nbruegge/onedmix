module onedmix_vmix_mytke
  use onedmix_variables
  use onedmix_eos
  use onedmix_utils
  use onedmix_io
  implicit none

  ! namelist parameters
  real*8 :: &
    c_k, c_eps, alpha_tke, mxl_min, kappaM_min, kappaM_max, cd, tke_surf_min, tke_min
  integer  :: &
    tke_mxl_choice
  logical :: &
    l_tke_active, only_tke, use_ubound_dirichlet, use_lbound_dirichlet

  ! TKE diagnostics
  real*8, dimension(:), allocatable :: &
                                       tke,          &
                                       Gimp_tke,      &
                                       Gexp_tke,      &
                                       tke_Tdif,     &
                                       tke_Tdis,     &
                                       tke_Tspr,     &
                                       tke_Tbpr,     &
                                       tke_Twin,     &
                                       tke_Ttot
  
  contains
!-------------------------------------------------------------------------------- 
  subroutine setup_vmix_mytke
    character(len=128)    :: fname
    namelist /tke_paras/                                                   &
        c_k, c_eps, alpha_tke, mxl_min, kappaM_min, kappaM_max, tke_mxl_choice   &
      , cd, tke_surf_min, tke_min, l_tke_active                                  &
      , only_tke, use_ubound_dirichlet, use_lbound_dirichlet
    ! read namelist or take standard parameters
    open(fid, file="./onedmix.nl", status="old", action='read')
    read(fid, nml=tke_paras)
    close(fid)

    allocate( tke(1:nz+1) ); tke=0.0
    allocate( Gimp_tke(1:nz+1) ); Gimp_tke=0.0
    allocate( Gexp_tke(1:nz+1) ); Gexp_tke=0.0

    allocate( tke_Tdif(1:nz+1) ); tke_Tdif=0.0
    allocate( tke_Tdis(1:nz+1) ); tke_Tdis=0.0
    allocate( tke_Tspr(1:nz+1) ); tke_Tspr=0.0
    allocate( tke_Tbpr(1:nz+1) ); tke_Tbpr=0.0
    allocate( tke_Twin(1:nz+1) ); tke_Twin=0.0
    allocate( tke_Ttot(1:nz+1) ); tke_Ttot=0.0

    ! read initial tke
    fname = trim(path_data) // "etke0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) tke
    close(fid)

  end subroutine setup_vmix_mytke

!-------------------------------------------------------------------------------- 
  subroutine calc_vmix_mytke
    real*8, dimension(nz+1) :: Gimp_tke_new, Gexp_tke_new, tke_exp
    real*8, dimension(nz+1) :: adiff, bdiff, cdiff, ldiag, mdiag, udiag
    real*8, dimension(nz+1) :: ke, Pr, cb
    real*8, dimension(nz+1) :: delta, Lmix, sqrttke
    real*8, dimension(nz+1) :: K_diss_v, P_diss_v
    real*8, dimension(nz+1) :: old_tke
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
    tke_exp = 0.0

    ! FIXME: Derive this!
    forc_rho_surf = 0.0
    !forc_tke_surf = 2e-5
    !taux = 2.4683998517220244E-004
    !tauy = 1.2652587216272979E-004
    !taux_act = 1.0550815445478266E-004
    !tauy_act = 9.9133200202586159E-006
    old_tke = tke

    sqrttke = sqrt(max(0d0,tke))
    Lmix = sqrt(2.0)*sqrttke/sqrt(max(1d-12, N2))
    ! !!! Do mixing length restriction stuff
    ! constrain mixing length scale as in MITgcm/OPA code (from pyOM)
    do k=2,nz+1
      Lmix(k) = MIN(Lmix(k), Lmix(k-1)+dzt(k-1) )
    end do 
    Lmix(1) = MIN( Lmix(1), mxl_min+dzt(1))
    do k=nz,1,-1
      Lmix(k) = MIN(Lmix(k), Lmix(k+1)+dzt(k))
    end do
    Lmix = max(Lmix, mxl_min)

    ! Pr=6.6*Ri with 1<Pr<10
    Pr = max(1.0, min(10.0,6.6*Ri) )
    Av = 0.
    kv = 0.
    !Av = Av + c_k * sqrt(tke) * Lmix
    !kv = kv + cb * sqrt(tke) * Lmix
    Av = min(kappaM_max, c_k*sqrttke*Lmix )
    kv = Av / Pr
    !Av = Av + Avb
    !kv = kv + kvb

    ! diffusion of tke
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
      Gimp_tke_new(k) = adiff(k)*tke(k-1) - bdiff(k)*tke(k) + cdiff(k)*tke(k+1)
    end do
    Gimp_tke_new(1)    = cdiff(1)*tke(2)     - bdiff(1)*tke(1)
    Gimp_tke_new(nz+1) = adiff(nz+1)*tke(nz) - bdiff(nz+1)*tke(nz+1)
    tke_Tdif = Gimp_tke_new
    ! add tke dissipation
    Gimp_tke_new = Gimp_tke_new - c_eps*sqrttke**3/Lmix
    tke_Tdis = -c_eps*sqrttke**3/Lmix

    ! explicit part
    ! add forcing
    P_diss_v   = kv*N2
    K_diss_v   = Av*S2
    P_diss_v(1) = -forc_rho_surf*grav/rho0
    tke_Tspr = K_diss_v
    tke_Tbpr = -P_diss_v
    Gexp_tke_new = K_diss_v-P_diss_v
    forc_tke_surf = cd * sqrt( (taux_act**2 + tauy_act**2)**(3./2.)  )
    Gexp_tke_new(1) = Gexp_tke_new(1) + forc_tke_surf/dzt(1)
    tke_Twin(1) = forc_tke_surf/dzt(1)

    ! Adams-Bashforth time stepping
    !write(*,*) 'tke_exp = ', tke_exp
    !write(*,*) 'tke = ', tke
    !write(*,*) 'Gimp_tke_new = ', Gimp_tke_new
    !write(*,*) 'Gimp_tke = ', Gimp_tke
    !write(*,*) 'Gexp_tke_new = ', Gexp_tke_new
    !write(*,*) 'Gexp_tke = ', Gexp_tke
    tke_exp = tke+dt*(1.d0-dimp)*((1.5+epsab)*Gimp_tke_new-(0.5+epsab)*Gimp_tke) &
                   +dt*            ((1.5+epsab)*Gexp_tke_new-(0.5+epsab)*Gexp_tke)
    Gimp_tke = Gimp_tke_new
    Gexp_tke = Gexp_tke_new

    ! implicit part
    ! FIXME: Add dissipation and forcing to mdiag
    udiag = -dt*dimp*cdiff
    mdiag = (1+dt*dimp*(bdiff + c_eps*sqrttke/Lmix))
    ldiag = -dt*dimp*adiff
    call solve_tridiag(ldiag, mdiag, udiag, tke_exp, tke, nz+1)
 
    ! set bounding values of tke
    if (tke(1) < 0.0) then
      tke(1) = 0.0
    end if
    tke(2:nz+1) = max(tke(2:nz+1), tke_min)
    tke_Ttot = (tke-old_tke)/dt

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
      write(*,*) 'sqrttke = ', sqrttke
      write(*,*) 'Av = ', Av
      write(*,*) 'delta*dt = ', delta*dt
      write(*,*) 'Gexp_tke_new = ', Gexp_tke_new
      !write(*,*) 'forc_tke_surf = ', forc_tke_surf
      write(*,*) 'ldiag = ', ldiag
      write(*,*) 'mdiag = ', mdiag
      write(*,*) 'udiag = ', udiag
      write(*,*) 'tke_exp = ', tke_exp
      write(*,*) 'tke = ', tke
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
    call save_variable(fprfx, tke, 'tke', 'turbulent kinetic energy', &
                       iostep, nz+1, 'm^2 / s^2')

    call save_variable(fprfx, tke_Tdif, 'tke_Tdif', &
                       'tke tendency by diffusion', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, tke_Tdis, 'tke_Tdis', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, tke_Tspr, 'tke_Tspr', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, tke_Tbpr, 'tke_Tbpr', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, tke_Twin, 'tke_Twin', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')

    call save_variable(fprfx, tke_Ttot, 'tke_Ttot', &
                       'tke tendency by dissipation', &
                       iostep, nz+1, 'm^2 / s^3')
  end subroutine write_snap_mytke

end module onedmix_vmix_mytke
