module onedmix_timeloop
  use onedmix_variables
  use onedmix_io
  use onedmix_eos
  use onedmix_vmix_mypp
  use onedmix_vmix_mytke
  use onedmix_vmix_myconst
  use onedmix_cvmix_tke
  implicit none
  contains
!-------------------------------------------------------------------------------- 
  subroutine timeloop
    integer :: k, l, ll
    real*8, dimension(nz) :: Gtemp, Gsalt, Guvel, Gvvel
    real*8, dimension(nz) :: Fexp_ptra, Fexp_temp, Fexp_salt, Fexp_uvel, Fexp_vvel
    real*8, dimension(nz) :: Fexp

    write(*,*) "================================================================================"
    write(*,*) "Start calculation..."
    write(*,*) "(nz, nt) = ", nz, nt
    write(*,*) "dt = ", dt
    write(*,*) "================================================================================"
    
    tact = 0
    do l=1,ntt
      do ll=1,nt

        !! explicit Adams-Bashforth time stepping
        !call calc_Gimp(temp, kv, dzw, dzt, nz, Gtemp)
        !temp = temp + ( (1.5+epsab)*Gtemp - (0.5+epsab)*Gtemp_old ) * dt
        !Gtemp_old = Gtemp
        !Gtemp = 0.0
  
        !call calc_Gimp(salt, kv, dzw, dzt, nz, Gsalt)
        !salt = salt + ( (1.5+epsab)*Gsalt - (0.5+epsab)*Gsalt_old ) * dt
        !Gsalt_old = Gsalt
        !Gsalt = 0.0
  
        !call calc_Gimp(uvel, Av, dzw, dzt, nz, Guvel)
        !uvel = uvel + ( (1.5+epsab)*Guvel - (0.5+epsab)*Guvel_old ) * dt
        !Guvel_old = Guvel
        !Guvel = 0.0
  
        !call calc_Gimp(vvel, Av, dzw, dzt, nz, Gvvel)
        !vvel = vvel + ( (1.5+epsab)*Gvvel - (0.5+epsab)*Gvvel_old ) * dt
        !Gvvel_old = Gvvel
        !Gvvel = 0.0
  
        Fexp_temp = 0.0
        Fexp_salt = 0.0
        Fexp_uvel = 0.0
        Fexp_vvel = 0.0
  
        ! --- interpolate surface forcing to current time
        ! (onedmix_timeloop/interp_forcing)
        call interp_forcing(q0, q0_act)
        Fexp_temp(1) = Fexp_temp(1) + q0_act/(cp*rho0)/dzw(1)
        call interp_forcing(emp, emp_act)
        Fexp_salt(1) = Fexp_salt(1) + emp_act/dzw(1)
        call interp_forcing(taux, taux_act)
        Fexp_uvel(1) = Fexp_uvel(1) + taux_act/dzw(1)
        call interp_forcing(tauy, tauy_act)
        Fexp_vvel(1) = Fexp_vvel(1) + tauy_act/dzw(1)

        !write(*,*) 'q0_act = ', q0_act
        !write(*,*) 'emp_act = ', emp_act
        !write(*,*) 'taux_act = ', taux_act
        !write(*,*) 'tauy_act = ', tauy_act

        !write(*,*) 'q0 = ', q0
        !write(*,*) 'emp = ', emp
        !write(*,*) 'taux = ', taux
        !write(*,*) 'tauy = ', tauy

        !if (ll==10) then
        !  stop
        !end if

        ! --- add bottom friction
        Fexp_uvel(nz) = Fexp_uvel(nz) &
          - bottomDragQuadratic*sqrt(0.5*(uvel(nz)**2+vvel(nz)**2))*uvel(nz)/dzw(nz)
        Fexp_vvel(nz) = Fexp_vvel(nz) &
          - bottomDragQuadratic*sqrt(0.5*(uvel(nz)**2+vvel(nz)**2))*vvel(nz)/dzw(nz)

        ! --- Coriolis force
        Fexp_uvel = Fexp_uvel + fCor*vvel
        Fexp_vvel = Fexp_vvel - fCor*uvel
  
        ! --- derive updated variable from explicite and implicite parts
        ! semi-implicit time stepping
        if (.true.) then
          ! (onedmix_timeloop/calc_Gimp2)
          call calc_Gimp2(temp, kv, dzw, dzt, nz, dt, dimp, epsab, temp, &
                            Gtemp_old, Gtemp_exp,                         &
                            Fexp_temp, 0)
          call calc_Gimp2(salt, kv, dzw, dzt, nz, dt, dimp, epsab, salt, &
                            Gsalt_old, Gsalt_exp,                         &
                            Fexp_salt, 0)
          call calc_Gimp2(uvel, Av, dzw, dzt, nz, dt, dimp, epsab, uvel, &
                            Guvel_old, Guvel_exp,                         &
                            Fexp_uvel, 0)
          call calc_Gimp2(vvel, Av, dzw, dzt, nz, dt, dimp, epsab, vvel, &
                            Gvvel_old, Gvvel_exp,                         &
                            Fexp_vvel, 0)
        end if
        Fexp_ptra = 0.0
        ! constant surface flux
        Fexp_ptra(1) = 1.0/dzw(1)
        ! restoring
        !Fexp_ptra(1) = 1.0/86400.0 * (1.0 - ptra(1))
        ! (onedmix_timeloop/calc_Gimp2)
        call calc_Gimp2(ptra, 1e-3*one_vec, dzw, dzt, nz, dt,  &
                         dimp, epsab, ptra,                     &
                         Gptra_old, Gptra_exp,                  &
                         Fexp_ptra, 0)
  
        ! --- derive actual density
        ! (onedmix_eos/calc_dens)
        call calc_dens(temp, salt, 0.d0, dens, nz)

        ! --- 
        ! (onedmix_timeloop/calc_vertical_gradients)
        call calc_vertical_gradients()
  
        ! --- derive mixing coefficients kv and Av
        if (mixing_scheme == 1) then
          ! (onedmix_vmix_mypp/calc_mypp)
          call calc_vmix_mypp() 
        elseif (mixing_scheme == 2) then
          ! (onedmix_vmix_mytke/calc_mytke)
          call calc_vmix_mytke()
        elseif (mixing_scheme == 3) then
          ! (onedmix_cvmix_tke/calc_cvmix_tke)
          call calc_cvmix_tke()
        elseif (mixing_scheme == 4) then
          ! (onedmix_vmix_myconst/calc_vmix_myconst)
          call calc_vmix_myconst()
        end if

        tact = tact+dt
      end do ! ll=1,nt
  
      ! --- model snapshot
      iostep = iostep + 1
      ! (onedmix_io/write_snapshot)
      call write_snapshot()

      if (mixing_scheme == 1) then
        ! (onedmix_vmix_mypp/write_snap_mypp)
        call write_snap_mypp()
      elseif (mixing_scheme == 2) then
        ! (onedmix_vmix_mypp/write_snap_mytke)
        call write_snap_mytke()
      elseif (mixing_scheme == 3) then
        ! (onedmix_cvmix_tke/write_snap_cvmix_tke)
        call write_snap_cvmix_tke()
      elseif (mixing_scheme == 4) then
        ! (onedmix_vmix_myconst/write_snap_myconst)
        call write_snap_myconst()
      end if
    end do ! l=1,ntt
  
  end subroutine timeloop

!-------------------------------------------------------------------------------- 
  subroutine calc_Gimp2(phi, kdiff, dzw, dzt, nz, dt, &
                         dimp, epsab, phi_new,         &
                         Gimp, Gexp,          &
                         Fexp, flag)
    implicit none
    
    integer, intent(in) :: nz
    real*8, intent(in), dimension(nz) :: phi, kdiff
    real*8, intent(in), dimension(nz) :: dzt, dzw
    real*8, intent(in), dimension(nz) :: Fexp
    real*8, intent(in) :: dimp, epsab, dt
    integer, intent(in) :: flag

    real*8, intent(inout), dimension(nz) :: Gimp, Gexp
    real*8, intent(out), dimension(nz) :: phi_new  

    real*8, dimension(nz) :: phi_exp, Gimp_new, Gexp_new 
    real*8, dimension(nz) :: adiff, bdiff, cdiff
    real*8, dimension(nz) :: ldiag, mdiag, udiag
    integer :: k
    ! time-stepping options
    ! Euler-Forward:      dimp=0.; epsab=-0.5
    ! Euler-Backward:     dimp=1.; epsab=-0.5
    ! Crank-Nicolson:     dimp=0.5; epsab=-0.5
    ! Adams-Bashforth:    dimp=0.; epsab=0.01
    ! mixed:              dimp=0.5; epsab=0.01

    adiff = 0.0
    bdiff = 0.0
    cdiff = 0.0
    ldiag = 0.0
    mdiag = 0.0
    udiag = 0.0
    Gimp_new = 0.0
    Gexp_new = 0.0
    phi_exp = 0.0

    do k=2,nz
      adiff(k) = kdiff(k)/(dzw(k)*dzt(k))
    end do
    adiff(1) = 0.0 ! adiff(1) does not enter diffusion matrix
    do k=2,nz-1
      bdiff(k) = kdiff(k)/(dzw(k)*dzt(k)) + kdiff(k+1)/(dzw(k)*dzt(k+1))
    end do
    bdiff(1)  = kdiff(2)/(dzw(1)*dzt(2)) ! dphi/dz(1)=0.
    bdiff(nz) = kdiff(nz)/(dzw(nz)*dzt(nz))  ! dphi/dz(nz+1)=0.
    do k=1,nz-1
      cdiff(k) = kdiff(k+1)/(dzw(k)*dzt(k+1))
    end do
    cdiff(nz) = 0.0 ! cdiff(nz) does not enter diffusion matrix

    ! implicit part
    do k=2,nz-1
      Gimp_new(k) = adiff(k)*phi(k-1) - bdiff(k)*phi(k) + cdiff(k)*phi(k+1)
    end do
    Gimp_new(1)  = cdiff(1)*phi(2) - bdiff(1)*phi(1)
    Gimp_new(nz) = adiff(nz)*phi(nz-1) - bdiff(nz)*phi(nz)

    !Gimp_new = 0.0
    !do k=2,nz-1
    !  Gimp_new(k) =       kdiff(k)/(dzw(k)*dzt(k))       * phi(k-1) &
    !                - (   kdiff(k)/(dzw(k)*dzt(k))                  &
    !                    + kdiff(k+1)/(dzw(k)*dzt(k+1)) ) * phi(k)   &
    !                +     kdiff(k+1)/(dzw(k)*dzt(k+1))   * phi(k+1) 
    !enddo
    !k = 1
    !Gimp_new(k) = - kdiff(k+1)/(dzw(k)*dzt(k+1))   * phi(k)   &
    !              + kdiff(k+1)/(dzw(k)*dzt(k+1))   * phi(k+1)
    !k = nz
    !Gimp_new(k) =   kdiff(k)/(dzw(k)*dzt(k))       * phi(k-1) &
    !              - kdiff(k)/(dzw(k)*dzt(k))       * phi(k)

    ! explicit part
    ! add forcing
    if (flag==1) then
      write(*,*) 'Fexp = ', Fexp
      write(*,*) 'phi = ', phi
    end if 
    Gexp_new  = Fexp

    ! Adams-Bashforth time stepping
    phi_exp = phi + dt*(1.d0-dimp)*( (1.5+epsab)*Gimp_new - (0.5+epsab)*Gimp ) &
                  + dt*            ( (1.5+epsab)*Gexp_new - (0.5+epsab)*Gexp )
    !Gimp_old = Gimp
    Gimp = Gimp_new
    Gexp = Gexp_new
    !phi_exp = phi + dt*(1.0-dimp)*Gimp

    ! implicit part
    udiag = -dt*dimp*cdiff
    mdiag = (1+dt*dimp*bdiff)
    ldiag = -dt*dimp*adiff

    if (.false.) then
      write(*,*) 'udiag = ', udiag(1:5)
      write(*,*) 'mdiag = ', mdiag(1:5)
      write(*,*) 'ldiag = ', ldiag(1:5)
    end if 

    call solve_tridiag(ldiag, mdiag, udiag, phi_exp, phi_new, nz)
    !phi_new = phi + dt*( (1.5+epsab)*Gimp_new - (0.5+epsab)*Gimp )
    !write(*,*) 'phi     = ', phi(1:5)
    !write(*,*) 'phi_exp = ', phi_exp(1:5)
    !write(*,*) 'phi_new = ', phi_new(1:5)

  end subroutine calc_Gimp2

!-------------------------------------------------------------------------------- 
  subroutine interp_forcing(forc, forc_interp)
    use onedmix_variables
    real*8, intent(in), dimension(nforc) :: forc
    real*8, intent(out) :: forc_interp
    integer :: nf
  
    ! linear interpolation
    nf = ceiling(tact/force_freq)
    !tl   = nf*force_freq
    !tlp1 = (nf+1)*force_freq
    !write(*,*) 'nf = ', nf
    !write(*,*) 'tact = ', tact
    !write(*,*) '---'
    !write(*,*) 'forc(nf+1) = ', forc(nf+1)
    !write(*,*) 'forc(nf)   = ', forc(nf)
    forc_interp = (forc(nf+1)-forc(nf))/force_freq * (tact - nf*force_freq) + forc(nf)
  end subroutine interp_forcing

!-------------------------------------------------------------------------------- 
  subroutine calc_vertical_gradients()
    real*8, dimension(nz+1) :: pint
    real*8                  :: dens_km1, dens_k
    integer :: k
    ! --- initialize/reset values
    N2=0.0
    S2=0.0
    uz=0.0
    vz=0.0
    Ri=0.0
    ! as in MPIOM beleg.f90 (l. 79)
    pint = 0.0001 * rho0 * (-zu)
    do k=2,nz
      !pint = 0.5*(pbcl(k-1)+pbcl(k))
      call potrho(temp(k-1), salt(k-1), pint(k), dens_km1) 
      call potrho(temp(k),   salt(k),   pint(k), dens_k) 
      N2(k) = -grav/rho0*(dens_km1-dens_k)/dzt(k)
      uz(k) = (uvel(k-1)-uvel(k))/dzt(k)
      vz(k) = (vvel(k-1)-vvel(k))/dzt(k)
      !write(*,*) 'dens_km1 = ', dens_km1
      !write(*,*) 'dens_k   = ', dens_k
    end do
    !write(*,*) 'dzt = ', dzt
    !write(*,*) 'N2 = ', N2
    !stop
    N2(1) = 0.0
    uz(1) = 0.0
    vz(1) = 0.0
    N2(nz+1) = 0.0
    uz(nz+1) = 0.0
    vz(nz+1) = 0.0
    S2 = uz**2+vz**2
    Ri = N2/max(S2,1e-12)
  end subroutine calc_vertical_gradients


!  subroutine calc_vert_grid()
!    use onedmix_variables
!    !real*8, dimension(nz), intent(in) :: dz
!    integer :: k
!    !dzw = dz
!    zu(1) = 0.0
!    do k=2,nz
!      zu(k) = zu(k-1) + dzw(k)
!    end do
!
!    do k=1,nz
!      zt(k) = zu(k)+0.5*dzw(k)
!    end do
!
!    dzt(1) = zt(1)
!    do k=2,nz
!      dzt(k) = zt(k)-zt(k-1)
!    end do
!
!    write(*,*) 'dzt = ', dzt
!  end subroutine calc_vert_grid


!-------------------------------------------------------------------------------- 
  ! FIXME: Is this needed anywhere?
  subroutine calc_Gdiff(phi, kdiff, dzw, dzt, nz, Gdiff)
  !subroutine calc_Gdiff(phi, kdiff, Gdiff)
  !  use onedmix_variables
  !  implicit none
    integer, intent(in) :: nz
    real*8, intent(in), dimension(nz) :: phi, kdiff
    real*8, intent(in), dimension(nz) :: dzt, dzw
    
    real*8, intent(out), dimension(nz) ::  Gdiff
    
    real*8, dimension(nz) :: phi_z, Fz
    integer :: k
    
    Gdiff = 0.0
    phi_z = 0.0
    Fz = 0.0
    
    do k=2,nz
      phi_z(k) = (phi(k-1)-phi(k))/dzt(k)
    end do
    
    Fz = kdiff*phi_z
    
    do k=1,nz-1
      Gdiff(k) = (Fz(k)-Fz(k+1))/dzw(k)
    end do
    Gdiff(nz) = Fz(nz)/dzw(nz)
  end subroutine calc_Gdiff
end module onedmix_timeloop
