module onedmix_setup
!
! 
!
  use onedmix_variables
  use onedmix_eos
  use onedmix_io
  use onedmix_vmix_mypp
  use onedmix_vmix_mytke 
  implicit none
contains

! -------------------------------------------------------------------------------- 
  subroutine setup_onedmix()
    ! --- read namelist "main"
    namelist/main/nz, nt, ntt, dt, dimp, epsab, Avb, kvb, rho0, grav, &
                  cal_type, cal_units, cal_origin, force_freq, nforc, &
                  mixing_scheme, fCor

    ! --- set default values
    ! FIXME: Define default values for every variable
    fCor = 0.0*1e-4

    ! --- read namelist
    if (.true.) then
      open(fid, file="./onedmix.nl", status="old", action='read')
      read(fid, nml=main)
      close(fid)
    end if
    ! check if enough forcing data is there
    if (nforc*force_freq < ntt*nt*dt) then
      write(*,*) '::: Error: Forcing not long enough for simulation length! :::'
      stop 2
    end if

    ! --- allocate basic variables 
    ! (onedmix_setup/allocate_vars)
    call allocate_vars() 

    ! --- read initial data and forcing
    ! (onedmix_setup/read_input_data)
    call read_input_data()

    ! --- calculate initial density
    ! (onedmix_eos/calc_dens)
    call calc_dens(temp, salt, 0.d0, dens, nz)

    ! --- write initial output file 
    ! (onedmix_io/write_snapshot)
    call write_snapshot()

    ! --- setup mixing scheme
    if (mixing_scheme == 1) then
      ! (onedmix_vmix_mypp/setup_mypp)
      call setup_vmix_mypp() 
      call write_snap_mypp()
    elseif (mixing_scheme == 2) then
      ! (onedmix_vmix_mytke/setup_mytke)
      call setup_vmix_mytke()
      call write_snap_mytke()
    end if

  end subroutine setup_onedmix
  
! -------------------------------------------------------------------------------- 
  subroutine allocate_vars()
    !use variables

    allocate( dzw(1:nz) ); dzw=0.0
    allocate( dzt(1:nz+1) ); dzt=0.0
    allocate( zt(1:nz) ); zt=0.0
    allocate( zu(1:nz+1) ); zu=0.0

    allocate( temp(1:nz) ); temp=0.0
    allocate( salt(1:nz) ); salt=0.0
    allocate( dens(1:nz) ); dens=0.0
    allocate( pbcl(1:nz) ); pbcl=0.0
    allocate( uvel(1:nz) ); uvel=0.0
    allocate( vvel(1:nz) ); vvel=0.0
    allocate( ptra(1:nz) ); ptra=0.0
    allocate( Gtemp_old(1:nz) ); Gtemp_old=0.0
    allocate( Gsalt_old(1:nz) ); Gsalt_old=0.0
    allocate( Guvel_old(1:nz) ); Guvel_old=0.0
    allocate( Gvvel_old(1:nz) ); Gvvel_old=0.0
    allocate( Gptra_old(1:nz) ); Gptra_old=0.0
    allocate( Gtemp_exp(1:nz) ); Gtemp_exp=0.0
    allocate( Gsalt_exp(1:nz) ); Gsalt_exp=0.0
    allocate( Guvel_exp(1:nz) ); Guvel_exp=0.0
    allocate( Gvvel_exp(1:nz) ); Gvvel_exp=0.0
    allocate( Gptra_exp(1:nz) ); Gptra_exp=0.0

    allocate( kv(1:nz+1) ); kv=0.0
    allocate( Av(1:nz+1) ); Av=0.0

    allocate( N2(1:nz+1) ); N2=0.0
    allocate( S2(1:nz+1) ); S2=0.0
    allocate( uz(1:nz+1) ); uz=0.0
    allocate( vz(1:nz+1) ); vz=0.0
    allocate( Ri(1:nz+1) ); Ri=0.0

    allocate( zero_vec(1:nz) ); zero_vec=0.0
    allocate( one_vec(1:nz) ); one_vec=1.0

    allocate( forc_time(1:nforc) ); forc_time=0.0
    allocate( q0(1:nforc) ); q0=0.0
    allocate( emp(1:nforc) ); emp=0.0
    allocate( taux(1:nforc) ); taux=0.0
    allocate( tauy(1:nforc) ); tauy=0.0

    lallocate = .true. 
  end subroutine allocate_vars

! -------------------------------------------------------------------------------- 
  subroutine read_input_data()
    !use variables
    character(len=128)    :: fname
    fname = trim(path_data) // "zt.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) zt
    close(fid)
    fname = trim(path_data) // "zu.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) zu
    close(fid)
    fname = trim(path_data) // "dzt.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) dzt
    close(fid)
    fname = trim(path_data) // "dzw.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) dzw
    close(fid)
    fname = trim(path_data) // "u0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) uvel
    close(fid)
    fname = trim(path_data) // "v0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) vvel
    close(fid)
    fname = trim(path_data) // "t0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) temp
    close(fid)
    fname = trim(path_data) // "s0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) salt
    close(fid)
    fname = trim(path_data) // "ptr0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) ptra
    close(fid)
    fname = trim(path_data) // "etke0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) Etke
    close(fid)
    fname = trim(path_data) // "forc_time.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) forc_time
    close(fid)
    fname = trim(path_data) // "forc_q0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) q0
    close(fid)
    fname = trim(path_data) // "forc_emp.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) emp
    close(fid)
    fname = trim(path_data) // "forc_taux.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) taux
    close(fid)
    fname = trim(path_data) // "forc_tauy.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) tauy
    close(fid)
  end subroutine

! -------------------------------------------------------------------------------- 
end module onedmix_setup

