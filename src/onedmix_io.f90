module onedmix_io
!
! This module contains function to do model io
!
  use onedmix_variables

  implicit none
contains

!-------------------------------------------------------------------------------- 
  subroutine write_snapshot()
    character(len=20)     :: fprfx

    write(*,*) "Write iostep: ", iostep

    ! --- nz variables
    fprfx = 'onedmix_state       '
    !call save_variable( fprfx, uvel, iostep, nz,                    &
    !                    'uvel                ',                     &
    !                    'm / sec             ',                     &
    !                    'zonal velocity                          ', &
    !                  )! iostep, nz, 'm / sec')
    call save_variable(fprfx, uvel, 'uvel', 'zonal velocity', &
                       iostep, nz, 'm / sec')

    call save_variable(fprfx, vvel, 'vvel', 'meridional velocity', &
                       iostep, nz, 'm / sec')

    call save_variable(fprfx, temp, 'temp', 'temperature', &
                       iostep, nz, 'deg C')

    call save_variable(fprfx, salt, 'salt', 'salinity', &
                       iostep, nz, 'g / kg')

    call save_variable(fprfx, dens, 'dens', 'density', &
                       iostep, nz, 'kg / m^3')

    call save_variable(fprfx, ptra, 'ptra', 'passive tracer concentr.', &
                       iostep, nz, '')
    
    ! --- nz+1 variables
    call save_variable(fprfx, Av, 'Av', 'vertical viscosity', &
                       iostep, nz+1, 'm^2 / s')

    call save_variable(fprfx, kv, 'kv', 'vertical diffusivity', &
                       iostep, nz+1, 'm^2 / s')

    call save_variable(fprfx, N2, 'N2', 'buoyancy frequency', &
                       iostep, nz+1, '1 / s^2')

    call save_variable(fprfx, S2, 'S2', 'vert. velocity shear', &
                       iostep, nz+1, '1 / s^2')

    call save_variable(fprfx, Ri, 'Ri', 'Richardson number', &
                       iostep, nz+1, '')
    ! --- 1D variables
    ! FIXME: add units
    call save_variable(fprfx, (/q0_act/), 'q0_act', 'surface heat flux', &
                       iostep, 1, '')

    call save_variable(fprfx, (/emp_act/), 'emp_act', 'surface salt flux', &
                       iostep, 1, '')

    call save_variable(fprfx, (/taux_act/), 'taux_act', 'zonal wind stress', &
                       iostep, 1, '')

    call save_variable(fprfx, (/tauy_act/), 'tauy_act', 'merid. wind stress', &
                       iostep, 1, '')

  end subroutine write_snapshot

!-------------------------------------------------------------------------------- 

  subroutine save_variable(fprfx, var, varname, long_name, iostep, nrec, units) 
    character(len=*)     :: fprfx
    character(len=*)     :: varname
    character(len=*)     :: long_name
    integer              :: iostep
    integer              :: nrec
    character(len=*)     :: units

    integer               :: recnum
    integer               :: bytes  = 4
    character(len=24)     :: endian = "big_endian"
    character(len=10)     :: tstepstr
    character(len=128)    :: fname

    ! FIXME: Does it work to use (:) without specifying exact dimension?
    real*8, dimension(:)  :: var


    ! --- tstepstr
    write(tstepstr,"(I10.10)") iostep 

    ! --- write data-file
    !fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".data"
    fname = trim(path_data) // '/out/' // &
            trim(varname) // "." // tstepstr // ".data"
    !write(*,*) trim(path_data)
    !write(*,*) trim(varname)
    !write(*,*) tstepstr
    !stop
    !write(*,*) varname, tstepstr, nrec
    open( unit=fid, file=fname, form='unformatted', status='replace', &
          access='direct', recl=bytes*nrec, convert=endian )
    recnum = 1
    write(fid, rec=recnum) sngl(var)
    close(fid)

    ! --- write meta-file
    !fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".meta"
    fname = trim(path_data) // '/out/' // &
            trim(varname) // "." // tstepstr // ".meta"
    open( unit=fid, file=fname, status='replace' )
    write(fid, *) 'name         = ', varname
    write(fid, *) 'long_name    = ', long_name
    write(fid, *) 'iostep       = ', iostep
    write(fid, *) 'nrec         = ', nrec
    write(fid, *) 'units        = ', units
    close(fid)
  end subroutine save_variable

!!-------------------------------------------------------------------------------- 
!  subroutine write_snapshot_old()
!    integer               :: recnum
!    integer               :: bytes  = 4
!    character(len=24)     :: endian = "big_endian"
!    character(len=20)     :: fprfx
!    character(len=10)     :: tstepstr
!    character(len=128)    :: fname
!    real*8, dimension(nz) :: outdata
!
!    write(*,*) "Write iostep: ", iostep
!
!    fprfx = 'onedmix_2D'
!    fname = trim(path_data) // trim(fprfx) // "_varlist.txt"
!    open( unit=fid+1, file=fname, status='replace' )
!
!    write(tstepstr,"(I10.10)") iostep 
!    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".data"
!    open( unit=fid, file=fname, form='unformatted', status='replace', &
!          access='direct', recl=bytes*nz, convert=endian )
!
!    recnum =1
!
!    write(fid+1, *) 'uvel'
!    write(fid, rec=recnum) sngl(uvel)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'vvel'
!    write(fid, rec=recnum) sngl(vvel)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'temp'
!    write(fid, rec=recnum) sngl(temp)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'salt'
!    write(fid, rec=recnum) sngl(salt)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'dens'
!    write(fid, rec=recnum) sngl(dens)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'ptra'
!    write(fid, rec=recnum) sngl(ptra)
!    recnum = recnum + 1
!
!    close(fid)
!    close(fid+1)
!
!    ! all variables that have length nz+1
!    fprfx = 'onedmix_2Dp1'
!    fname = trim(path_data) // trim(fprfx) // "_varlist.txt"
!    open( unit=fid+1, file=fname, status='replace' )
!
!    write(tstepstr,"(I10.10)") iostep 
!    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".data"
!    open( unit=fid, file=fname, form='unformatted', status='replace', &
!          access='direct', recl=bytes*(nz+1), convert=endian )
!
!    recnum =1
!
!    write(fid+1, *) 'Av'
!    write(fid, rec=recnum) sngl(Av)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'kv'
!    write(fid, rec=recnum) sngl(kv)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'N2'
!    write(fid, rec=recnum) sngl(N2)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'S2'
!    write(fid, rec=recnum) sngl(S2)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'Ri'
!    write(fid, rec=recnum) sngl(Ri)
!    recnum = recnum + 1
!
!    ! FIXME: Put this somehow to 1dmix_vmix_mytke.f90
!    if (mixing_scheme==2) then
!      write(fid+1, *) 'Etke'
!      write(fid, rec=recnum) sngl(Etke)
!      recnum = recnum + 1
!
!      write(fid+1, *) 'TEtke_dif'
!      write(fid, rec=recnum) sngl(TEtke_dif)
!      recnum = recnum + 1
!
!      write(fid+1, *) 'TEtke_dis'
!      write(fid, rec=recnum) sngl(TEtke_dis)
!      recnum = recnum + 1
!
!      write(fid+1, *) 'TEtke_spr'
!      write(fid, rec=recnum) sngl(TEtke_spr)
!      recnum = recnum + 1
!
!      write(fid+1, *) 'TEtke_bpr'
!      write(fid, rec=recnum) sngl(TEtke_bpr)
!      recnum = recnum + 1
!
!      write(fid+1, *) 'TEtke_tau'
!      write(fid, rec=recnum) sngl(TEtke_tau)
!      recnum = recnum + 1
!
!      write(fid+1, *) 'TEtke_tot'
!      write(fid, rec=recnum) sngl(TEtke_tot)
!      recnum = recnum + 1
!    end if
!
!    close(fid)
!    close(fid+1)
!    
!    ! 1D variables (time-series of scalars)
!    fprfx = 'onedmix_1D'
!    fname = trim(path_data) // trim(fprfx) // "_varlist.txt"
!    open( unit=fid+1, file=fname, status='replace' )
!
!    write(tstepstr,"(I10.10)") iostep 
!    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".data"
!    open( unit=fid, file=fname, form='unformatted', status='replace', &
!          access='direct', recl=bytes, convert=endian )
!    recnum = 1
!
!    write(fid+1, *) 'q0_act'
!    write(fid, rec=recnum) sngl(q0_act)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'emp_act'
!    write(fid, rec=recnum) sngl(emp_act)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'taux_act'
!    write(fid, rec=recnum) sngl(taux_act)
!    recnum = recnum + 1
!
!    write(fid+1, *) 'tauy_act'
!    write(fid, rec=recnum) sngl(tauy_act)
!    recnum = recnum + 1
!
!    close(fid)
!    close(fid+1)
!  end subroutine write_snapshot_old

end module onedmix_io
