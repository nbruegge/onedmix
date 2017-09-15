
module cvmix_idemix
! This module contains the main computations of the IDEMIX 1 parameterization (described in "A Global Model for the Diapycnal
! Diffusivity Induced by Internal Gravity Waves", Olbers&Eden 2013) of Internal wave energy and its dissipation
! 
! @see
!
! @see
!
!
!! @author Hannah Kleppin, MPIMET/University of Hamburg
!! @author Oliver Gutjahr, MPIMET
!!
!! @par Copyright
!! 2002-2013 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! FIXME: lateral diffusion within module?

use cvmix_kinds_and_types,    only : cvmix_r8,                     &
                                      CVMIX_OVERWRITE_OLD_VAL,     &
                                      CVMIX_SUM_OLD_AND_NEW_VALS,  &
                                      CVMIX_MAX_OLD_AND_NEW_VALS,  &
                                      cvmix_data_type,             &
                                      cvmix_PI,                    & 
                                      cvmix_global_params_type

use cvmix_utils,              only : cvmix_update_tke, solve_tridiag


implicit none
private 
save


!public member functions

public :: init_idemix
public :: cvmix_coeffs_idemix
public :: gofx2  ! fixme: only used by IDEMIX public?
public :: hofx1  ! fixme: public?

!=================================================================================
!---------------------------------------------------------------------------------
! Interface to call the IDEMIX parameterization
!---------------------------------------------------------------------------------

interface cvmix_coeffs_idemix
    module procedure integrate_idemix  ! calculation ! FIXME: rename in cvmix_coeffs_low..
    !module procedure idemix_wrap       ! necessary to handle old/new values and to hand over user_defined constants
end interface cvmix_coeffs_idemix

!---------------------------------------------------------------------------------
! Interface to put values to IDEMIX variables
!---------------------------------------------------------------------------------

! FIXME: rename procedures
interface idemix_put
    module procedure vmix_tke_put_idemix_int
    module procedure vmix_tke_put_idemix_real
end interface idemix_put

!=================================================================================

! types for Idemix
type, public :: idemix_type
private


! through idemix_put this parameters are set to either default or user defined values
 real(cvmix_r8) ::  &
   tau_v            ,& ! time scale for vertical symmetrisation (unit?)
   !tau_h           ,& ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (unit?)
   gamma            ,& ! constant of order one derived from the shape of the spectrum in m space (unit?) (Olbers & Eden XXXX)
   jstar            ,& ! spectral bandwidth in modes (unit?)
   mu0                 ! dissipation parameter (unit?)

! Flag for how to update old values 
! Note: We don't need max or sum option
 integer :: handle_old_vals
end type idemix_type

type(idemix_type), target :: idemix_constants_saved 

contains

!=================================================================================

subroutine init_idemix(tau_v, gamma,jstar,mu0,handle_old_vals,idemix_userdef_constants)

! This subroutine sets user or default values for IDEMIX parameters

real(cvmix_r8),optional, intent(in) ::      &
  tau_v                                    ,& ! 
  gamma                                    ,& !
  jstar                                    ,& !
  mu0

type(idemix_type), intent(inout),target, optional :: idemix_userdef_constants

integer,intent(in),optional :: handle_old_vals

! FIXME: not sure about the allowed ranges for idemix parameters, default values confirm with pyOM testcases
if (present(tau_v)) then
  if(tau_v.lt.1.d0*86400.0 .or. tau_v .gt. 100.d0*86400.0) then
    print*, "ERROR:tau_v can only be allowed_range"
    stop 1
  end if
  call idemix_put('tau_v', tau_v, idemix_userdef_constants)
else
  call idemix_put('tau_v',1.d0*86400.0 , idemix_userdef_constants)
end if

! FIXME: only neccessary for lateral diffusion! 
! if (present(tau_h)) then
!   if(tau_h.lt. number .or. tau_h .gt. number) then
!     print*, "ERROR:tau_h can only be allowed_range"
!     stop 1
!   end if
!   call idemix_put('tau_h', tau_h, idemix_userdef_constants)
! else
!   call idemix_put('tau_h', 15.0*86400.0, idemix_userdef_constants)
! end if

if (present(gamma)) then
  if(gamma.lt. 1.d0 .or. gamma .gt. 3.d0) then
    print*, "ERROR:gamma can only be allowed_range"
    stop 1
  end if
  call idemix_put('gamma', gamma, idemix_userdef_constants)
else
  call idemix_put('gamma', 1.57d0, idemix_userdef_constants)
end if

if (present(jstar)) then
  if(jstar.lt. 5.d0 .or. jstar .gt. 15.d0) then
    print*, "ERROR:jstar can only be allowed_range"
    stop 1
  end if
  call idemix_put('jstar', jstar, idemix_userdef_constants)
else
  call idemix_put('jstar', 10.d0 , idemix_userdef_constants)
end if

if (present(mu0)) then
  if(mu0.lt. 1.d0 .or. mu0 .gt. 3.d0) then
    print*, "ERROR: mu0 can only be allowed_range"
    stop 1
  end if
  call idemix_put('mu0', mu0, idemix_userdef_constants)
else
  call idemix_put('mu0', 4.d0/3.0 , idemix_userdef_constants)
end if

if (present(handle_old_vals)) then
  if(handle_old_vals.lt. 1 .or. handle_old_vals.gt. 3 ) then
    print*, "ERROR:handle_old_vals can only be 1 to 3"
    stop 1
  end if
  call idemix_put('handle_old_vals', handle_old_vals, idemix_userdef_constants)
else
  call idemix_put('handle_old_vals', 1, idemix_userdef_constants)
end if

end subroutine init_idemix

!=================================================================================

subroutine idemix_wrap(Vmix_vars, idemix_userdef_constants)

! This subroutine is necessary to handle old/new values and to hand over the IDEMIX parameter set in previous subroutine
! This subroutine should be called from calling ocean model or driver

type(idemix_type), intent(in), optional, target :: idemix_userdef_constants

type(cvmix_data_type), intent(inout) :: Vmix_vars

! NEW values
real(cvmix_r8), dimension(Vmix_vars%nlev+1) ::    &
  new_E_iw                                       ,& !
  cvmix_int_1                                       ,& !
  cvmix_int_2                                       ,& !
  cvmix_int_3                                       ,& !
  iwe_Ttot                                          ,&
  iwe_Tdif                                          ,&
  iwe_Tdis                                          ,&
  iwe_Tsur                                          ,&
  iwe_Tbot                                          ,&
  new_KappaM                                        ,&
  new_KappaH                                        ,&
  c0                                                ,&
  new_iw_diss                                       ! 

integer ::                                        &
  nlev                                           ,& !
  max_nlev                                          !
! nils
integer :: i,j, tstep_count

type(idemix_type), pointer :: idemix_constants_in


idemix_constants_in => idemix_constants_saved
if (present(idemix_userdef_constants)) then
  idemix_constants_in => idemix_userdef_constants
end if

nlev = Vmix_vars%nlev
max_nlev = Vmix_vars%max_nlev


! call to actual computation of IDEMIX parameterization
! nils: added i,j
write(*,*) 'I am wrapping'
stop
call cvmix_coeffs_idemix( &
                         i = i, &
                         j = j, &
                         tstep_count = tstep_count, &
                         iw_diss_out     = new_iw_diss,               &
                         E_iw_out        = new_E_iw,                  &
                         KappaM_out      = new_KappaM,                &
                         KappaH_out      = new_KappaH,                &
                         cvmix_int_1     = cvmix_int_1,               &
                         cvmix_int_2     = cvmix_int_2,               &
                         cvmix_int_3     = cvmix_int_3,               &
                         iwe_Ttot        = iwe_Ttot,                  &
                         iwe_Tdif        = iwe_Tdif,                  &
                         iwe_Tdis        = iwe_Tdis,                  &
                         iwe_Tsur        = iwe_Tsur,                  &
                         iwe_Tbot        = iwe_Tbot,                  &
                         dzw             = Vmix_vars%dzw,             &
                         dzt             = Vmix_vars%dzt,             &
                         nlev            = nlev,                      &
                         max_nlev        = max_nlev,                  &
                         old_E_iw        = Vmix_vars%E_iw,            &
                         old_iw_diss     = Vmix_vars%iw_diss,         &
                         Nsqr            = Vmix_vars%Nsqr_iface,      &
                         Ssqr            = Vmix_vars%Ssqr_iface,      &
                         forc_iw_surface = Vmix_vars%forc_iw_surface, &
                         dtime           = Vmix_vars%dtime,           &
                         forc_iw_bottom  = Vmix_vars%forc_iw_bottom,  &
                         Kappa_GM        = Vmix_vars%Kappa_GM,        &
                         coriolis        = Vmix_vars%coriolis,        &
                         alpha_c         = Vmix_vars%alpha_c,         &
                         c0              = c0,                        &
                         idemix_userdef_constants = idemix_userdef_constants)

! update Vmix_vars to new values
! FIXME: nils: This should probably be cvmix_update_idemix. Hoever, it is not used
! anyway.
call cvmix_update_tke(idemix_constants_in%handle_old_vals,            &
                      nlev,                                           &
                      iw_diss_out = Vmix_vars%iw_diss,                &
                      new_iw_diss = new_iw_diss,                      &
                      E_iw_out    = Vmix_vars%E_iw,                   &
                      new_E_iw    = new_E_iw)

end subroutine idemix_wrap

!=================================================================================

subroutine integrate_idemix( &
                            i, j, &
                            tstep_count, &
                            iw_diss_out,                         & 
                            E_iw_out,                            &
                            KappaM_out,                          &
                            KappaH_out,                          &
                            cvmix_int_1,                         &
                            cvmix_int_2,                         &
                            cvmix_int_3,                         &
                            iwe_Ttot,                            &
                            iwe_Tdif,                            &
                            iwe_Tdis,                            &
                            iwe_Tsur,                            &
                            iwe_Tbot,                            &
                            dzw,                                 &
                            dzt,                                 &
                            nlev,                                &
                            max_nlev,                            &
                            old_E_iw,                            &
                            old_iw_diss,                         &
                            Nsqr,                                &
                            Ssqr,                                &
                            forc_iw_surface,                     &
                            dtime,                               &
                            forc_iw_bottom,                      &
                            Kappa_GM,                            & 
                            coriolis,                            &
                            alpha_c,                             &
                            c0,                                  &
                            idemix_userdef_constants)

! This subroutine contains the actual computation of IDEMIX
  
 ! Local varaibles
 type(idemix_type), intent(in), optional, target :: idemix_userdef_constants

 integer, intent(in)                                     ::      &
   nlev                                                         ,&
   max_nlev                                                         

 real(cvmix_r8), dimension(nlev+1), intent(inout)             :: &
    KappaM_out                                                  ,&
    KappaH_out

 ! nils
 integer, intent(in) :: i, j, tstep_count

 real(cvmix_r8), dimension(nlev+1), intent(in)           ::      &
   dzw

 real(cvmix_r8), dimension(nlev+1), intent(in)           ::      &
   Nsqr                                                         ,&
   old_E_iw                                                     ,&
   old_iw_diss                                                  ,& 
   dzt                                                             !

 ! diagnostics
 real(cvmix_r8), dimension(nlev+1), intent(out) ::               &
   iw_diss_out                                                  ,& 
   E_iw_out                                                     ,&
   cvmix_int_1                                                  ,&
   cvmix_int_2                                                  ,&
   cvmix_int_3                                                  ,&
   iwe_Ttot                                                     ,&
   iwe_Tdif                                                     ,&
   iwe_Tdis                                                     ,&
   iwe_Tsur                                                     ,&
   iwe_Tbot                                                     ,&
   c0                                                           ,&
   alpha_c

! not implemented currently, could be added once energy conserving linking
! between diff. parameterizations is desired
 real(cvmix_r8),dimension(nlev+1), intent(in), optional  ::      &
   Kappa_GM                                                     ,& ! 
   Ssqr

 real(cvmix_r8), intent(in)                              ::      & 
   forc_iw_bottom                                               ,& !
   !forc_iw_surface                                              ,& !
   dtime                                                        ,& !
   coriolis                                                        !
 real(cvmix_r8) ::  forc_iw_surface

 integer                                                 ::      &
   k, ks, ke, n

 ! input to the tri-diagonal solver
 real(cvmix_r8), dimension(nlev+1)                       ::      &
   a_dif                                                        ,& !
   b_dif                                                        ,& !
   c_dif                                                        ,& !
   a_tri                                                        ,& !
   b_tri                                                        ,& !
   c_tri                                                        ,& !
   d_tri                                                        ,& !
   delta                                                        ,& !
   !c0                                                           ,& ! mean vertical group velocity                      
   v0                                                           ,& ! mean lateral group velocity
   maxE_iw                                                      ,& ! 
   forc                                                            ! forcing, (currently not used)

 !IDEMIX parameters
 real(cvmix_r8)                                          ::      & 
   cstar                                                        ,& ! modal gravity wave speed
   fxa                                                          ,& !
   tau_v                                                        ,& !
!  tau_h                                                        ,& !
   gamma                                                        ,& !
   jstar                                                        ,& !
   mu0                                                          ,& !
   bN0                                                             !

 type(idemix_type), pointer ::idemix_constants_in

 ! initialize diagnostics
 iwe_Ttot = 0.0
 iwe_Tdif = 0.0
 iwe_Tdis = 0.0
 iwe_Tsur = 0.0
 iwe_Tbot = 0.0

 idemix_constants_in => idemix_constants_saved
 if (present(idemix_userdef_constants)) then
   idemix_constants_in => idemix_userdef_constants
 end if

 ! set idemix_constants locally
 tau_v = idemix_constants_in%tau_v
 gamma = idemix_constants_in%gamma
 mu0   = idemix_constants_in%mu0
 jstar = idemix_constants_in%jstar

 ! FIXME: nils: What should we do with the comment below?
 ! include "mass.include"  ! include this on AIX which does not know function acosh, also link with -lmass

 ! calculate cstar from OE13 Eq. (13)
 bN0=0.0
 do k=2,nlev
   bN0 = bN0 + max(0d0,Nsqr(k))**0.5*dzw(k) 
   !print*,"bn0=",bN0
 enddo
 ! FIXME: nils: check if Nsqr(1)!=0; should it be?
 bN0 = bN0 + max(0d0,Nsqr(1))**0.5*0.5*dzw(1) 
 cstar = max(1d-2,bN0/(cvmix_PI*jstar) )
    
 ! calculate vertical and horizontal representative group velocities c0 and v0
 ! c0: OE13 Eq. (13) and v0: OE13 Eq. (A9)
 ! alpha_c iwe**2: dissipation of internal wave energy (OE13 Eq. (15))
 do k=1,nlev+1
   ! FIXME: Carste: gofx2 bug 
   fxa = max(0d0,Nsqr(k))**0.5/(1d-22 + abs(coriolis) )
   ! print*,"fxa=",fxa
   c0(k)=max(0d0, gamma*cstar*gofx2(fxa) )
   v0(k)=max(0d0, gamma*cstar*hofx1(fxa))
   alpha_c(k) = max( 1d-4, mu0*acosh(max(1d0,fxa))*abs(coriolis)/cstar**2 )
 enddo

 ! FIXME: nils: remove next line
 !c0 = 1.0

!---------------------------------------------------------------------------------
! once energy conserving linking between the diff. parameterizations is required different forcings for IDEMIX could be added here
!---------------------------------------------------------------------------------
 
 forc(:)=0.d0
! ! all optional
! k_diss_gm=Kappa_GM*Ssqr
! SHORTCUT WITHOUT EKE MODEL:
! later on EKE dissipation could be added here
! forc= K_diss_h + k_diss_gm - P_diss_skew
! 
! ! K_diss_h is dissipation by harmonic/biharmonic friction
! if (.not. bottom_friction_tke) forc = forc + K_diss_bot
! 
! ! Maybe treatment of surf and bottom forcing could be done here as well

!---------------------------------------------------------------------------------

 ! prevent negative dissipation of IW energy
 ! FIXME: Carsten thinks we don't need this
 maxE_iw = max(0D0, old_E_iw)


 ! vertical diffusion and dissipation is solved implicitely 
 !---------------------------------------------------------------------------------
 ! assignment of tridiagonal matrix
 !---------------------------------------------------------------------------------
 ! |b1 c1 0  0  0  | (E1) = (d1)
 ! |a2 b2 c2 0  0  | (E2) = (d2)
 ! |0  a3 b3 c3 0  | (E3) = (d3)
 ! |0  0  a4 b4 c4 | (E4) = (d4)
 ! |0  0  0  an bn | (En) = (dn)
 !
 ! d1 = diss_1 + surf_forc 
 ! dn = diss_n + bott_forc 
 ! 
   
 ! vertical flux
 do k=1,nlev
  !delta(k) = dtime*tau_v/dzw(k-1)*0.5_cvmix_r8*(c0(k)+c0(k-1))
  delta(k) = tau_v/dzw(k) * 0.5*(c0(k)+c0(k+1))
 !if (i==45 .and. j==10) then
 !   write(*,*) 'k, nlev = ', k, nlev
 !   write(*,*) 'delta(k) = ', delta(k)
 !   write(*,*) 'dtime = ', dtime, tau_v
 !   write(*,*) 'dzw(k-1) ', dzw(k-1)
 !   write(*,*) 'c0 = ', c0(k), c0(k-1)
 !   write(*,*) 'dtime*... = ', dtime*tau_v/dzw(k-1)
 !   write(*,*) 0.5_cvmix_r8*(c0(k)+c0(k-1))
 ! end if
 enddo
 !delta(1)=0.0  
 delta(nlev+1) = 0.0          ! delta(nlev+1) is never used

 ! -- a -- 
 do k=2,nlev+1
   !a_tri(k) = - delta(k)*c0(k-1)/dzt(k)
   a_dif(k) = delta(k-1)*c0(k-1)/dzt(k)
 enddo
 a_dif(1) = 0.0 ! not part of the diffusion matrix, thus value is arbitrary

 ! -- b -- 
 do k=2,nlev
   !b_tri(k) = 1 + delta(k)*c0(k)/dzt(k) + delta(k+1)*c0(k)/dzt(k) + dtime*alpha_c(k)*maxE_iw(k)
   b_dif(k) = (delta(k-1)*c0(k)+delta(k)*c0(k))/dzt(k)
 enddo
 !k = 1
 !! FIXME: Check if really 0.5*dzt(1) should be taken!
 !b_tri(k) = 1 + delta(k+1)*c0(k)/(0.5*dzt(k)) + dtime*alpha_c(k)*maxE_iw(k)
 !k = nlev
 !b_tri(k) = 1 + delta(k)*c0(k)/dzt(k)         + dtime*alpha_c(k)*maxE_iw(k)

 ! Neumann boundary conditions
 k = 1
 b_dif(k) = delta(k)*c0(k)/dzt(k)
 k = nlev+1
 b_dif(k) = delta(k-1)*c0(k)/dzt(k)

 ! -- c-- 
 do k=1,nlev
   !c_tri(k) = - delta(k+1)/dzt(k)*c0(k+1)
   c_dif(k) = delta(k)*c0(k+1)/dzt(k)
 enddo
 c_dif(nlev+1) = 0.0 ! not part of the diffusion matrix, thus value is arbitrary

 !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
 a_tri = -dtime*a_dif
 b_tri = 1+dtime*b_dif
 ! FIXME: nils: Should dissipation also be in first and last layer?
 b_tri(2:nlev) = b_tri(2:nlev) + dtime*alpha_c(2:nlev)*maxE_iw(2:nlev)
 c_tri = -dtime*c_dif
  
 ! -- d -- 
 d_tri(1:nlev+1) = old_E_iw(1:nlev+1) + dtime*forc(1:nlev+1)
 d_tri(nlev+1)   = d_tri(nlev+1)      + dtime*forc_iw_bottom/dzt(nlev+1) 
 d_tri(1)        = d_tri(1)           + dtime*forc_iw_surface/dzt(1)

 ! solve the tri-diag matrix 
 call solve_tridiag(a_tri, b_tri, c_tri, d_tri, E_iw_out, nlev+1)

 ! --- diagnose implicite tendencies (only for diagnostics)
 ! vertical diffusion of E_iw
 do k=2,nlev
   iwe_Tdif(k) = a_dif(k)*E_iw_out(k-1) - b_dif(k)*E_iw_out(k) + c_dif(k)*E_iw_out(k+1)
 enddo
 k = 1
 iwe_Tdif(k) = - b_dif(k)*E_iw_out(k) + c_dif(k)*E_iw_out(k+1)
 k = nlev+1
 iwe_Tdif(k) = a_dif(k)*E_iw_out(k-1) - b_dif(k)*E_iw_out(k)

 ! dissipation of IW energy
 iw_diss_out = alpha_c * maxE_iw * E_iw_out
 iwe_Tdis = -alpha_c * maxE_iw * E_iw_out

 iwe_Tsur(1)      = forc_iw_surface/dzt(1) 
 iwe_Tbot(nlev+1) = forc_iw_bottom/dzt(nlev+1)

 iwe_Ttot = (E_iw_out-old_E_iw)/dtime

 ! derive diffusivity and viscosity by using Osbourne-Cox relation
 KappaH_out = 0.0
 KappaM_out = 0.0
 do k=2,nlev
   KappaH_out(k) =  0.2/(1.0+0.2) * iw_diss_out(k) / max(1d-12, Nsqr(k))
   KappaH_out(k) = max(1e-9, KappaH_out(k))
   KappaH_out(k) = min(1.0, KappaH_out(k))
   KappaM_out(k) =  10.0 * KappaH_out(k)
 enddo

 !---------------------------------------------------------------------------------
 cvmix_int_1 = Nsqr
 cvmix_int_2 = Nsqr 
 cvmix_int_3 = c0

 ! nils
 if (.false.) then
 if (i==45 .and. j==10) then
    write(*,*) ' ===================== '

    write(*,*) 'dtime = ', dtime
    write(*,*) 'delta = ', delta
    write(*,*) 'dzw = ', dzw
    write(*,*) 'c0 = ', c0
    write(*,*) 'a_tri = ', a_tri
    write(*,*) 'b_tri = ', b_tri
    write(*,*) 'c_tri = ', c_tri
    write(*,*) 'd_tri = ', d_tri
    write(*,*) 'forc_iw_surface = ', forc_iw_surface
    write(*,*) 'E_iw_out = ', E_iw_out

    write(*,*) 'iwe_Ttot = ', iwe_Ttot
    write(*,*) 'iwe_Tdif = ', iwe_Tdif
    write(*,*) 'iwe_Tdis = ', iwe_Tdis
    write(*,*) 'iwe_Tsur = ', iwe_Tsur
    write(*,*) 'iwe_Tbot = ', iwe_Tbot
    write(*,*) 'iwe_Tres = ', iwe_Ttot-(iwe_Tdif+iwe_Tdis+iwe_Tsur+iwe_Tbot)

   write(*,*) 'tau_v = ', tau_v
   write(*,*) 'gamma = ', gamma
   write(*,*) 'jstar = ', jstar
   write(*,*) 'mu0 = ', mu0

   !stop
 endif
 endif

end subroutine integrate_idemix

!=================================================================================


function gofx2(x)
!=======================================================================
! a function g(x)	!from pyOM 
!=======================================================================
 implicit none
 real(cvmix_r8) :: gofx2,x,c
 real*8, parameter :: pi = 3.14159265358979323846264338327950588
 x=max(3d0,x)
 c= 1.-(2./pi)*asin(1./x)
 gofx2 = 2/pi/c*0.9*x**(-2./3.)*(1-exp(-x/4.3))
end function gofx2

function hofx1(x)
!=======================================================================
! a function h(x) 	!from pyOM
!=======================================================================
 implicit none
 real(cvmix_r8) :: hofx1,x
 real*8, parameter :: pi = 3.14159265358979323846264338327950588
 hofx1 = (2./pi)/(1.-(2./pi)*asin(1./x)) * (x-1.)/(x+1.)
end function hofx1

!=================================================================================

subroutine vmix_tke_put_idemix_real(varname,val,idemix_userdef_constants)

! This subroutine puts real values to IDEMIX variables
!IN
  character(len=*),          intent(in) :: varname
  real(cvmix_r8),            intent(in) :: val
!OUT   
  type(idemix_type), intent(inout), target, optional:: idemix_userdef_constants
  type(idemix_type), pointer :: idemix_constants_out

  idemix_constants_out=>idemix_constants_saved
  if (present(idemix_userdef_constants)) then
    idemix_constants_out=> idemix_userdef_constants
  end if
  
select case(trim(varname))

    case('tau_v') 
      idemix_constants_out%tau_v= val
    case('jstar') 
      idemix_constants_out%jstar= val
    case('gamma') 
      idemix_constants_out%gamma = val
    case('mu0') 
      idemix_constants_out%mu0 = val
   
    case DEFAULT
      print*, "ERROR:", trim(varname), " not a valid choice"
      stop 1
end select

end subroutine vmix_tke_put_idemix_real

!=================================================================================

subroutine vmix_tke_put_idemix_int(varname,val,idemix_userdef_constants)

! This subroutine puts integer values to IDEMIX variables
!IN
  character(len=*),           intent(in) :: varname
  integer,                    intent(in) :: val
!OUT   
  type(idemix_type), intent(inout), target, optional:: idemix_userdef_constants
  type(idemix_type), pointer :: idemix_constants_out

  idemix_constants_out=>idemix_constants_saved
  if (present(idemix_userdef_constants)) then
    idemix_constants_out=> idemix_userdef_constants
  end if

  select case(trim(varname))

    case('handle_old_vals')
      idemix_constants_out%handle_old_vals=val
    
  end select
    
end subroutine vmix_tke_put_idemix_int

!=================================================================================

end module cvmix_idemix 
