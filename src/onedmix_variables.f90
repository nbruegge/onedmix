module onedmix_variables
!
! This module works as a container for all important model variables.
!

  implicit none

  !integer :: int1=10
  !real*8  :: real1=10.

  ! namelist parameters
  integer             :: nz, nt, ntt
  real*8              :: dt
  real*8              :: dimp, epsab
  ! real*8             :: kv_min, kv_max, Av_min, Av_max ! FIXME: Do we need those?
  real*8              :: kvB, Avb  ! FIXME: Put these to PP scheme
  real*8              :: fCor ! Coriolis parameter
  real*8              :: bottomDragQuadratic ! quadratic bottom-drag coefficient (dimensionless)
  real*8              :: rho0 != 1025. kg/m^3
  real*8              :: cp   != 4000. J/(kg K)
  real*8              :: grav != 9.81  m/s^2
  character(len=128)  :: cal_type, cal_units, cal_origin
  real*8              :: force_freq
  integer             :: nforc
  integer             :: mixing_scheme
  logical             :: no_slip_bottom
  integer             :: eos_type

  integer, parameter :: nn=9
  logical :: lallocate=.false.
  real*8, dimension(nn) :: test
  
  integer :: iostep=0
  integer :: tstep_count=0

  real*8, allocatable, dimension(:) :: dzt,           &
                                       dzw,           &
                                       zt,            &
                                       zu

  real*8, allocatable, dimension(:) :: temp,          &
                                       salt,          &
                                       dens,          &
                                       pbcl,          &
                                       uvel,          &
                                       vvel,          &
                                       ptra,          &
                                       Gtemp_old,     &
                                       Gsalt_old,     &
                                       Guvel_old,     &
                                       Gvvel_old,     &
                                       Gptra_old,     &
                                       Gtemp_exp,     &
                                       Gsalt_exp,     &
                                       Guvel_exp,     &
                                       Gvvel_exp,     &
                                       Gptra_exp

  real*8, dimension(:), allocatable :: &
                                       kv,            &
                                       Av
  real*8, allocatable, dimension(:) :: N2, S2, uz, vz, Ri 
  real*8, allocatable, dimension(:) :: zero_vec, one_vec
  real*8, allocatable, dimension(:) :: forc_time, q0, emp, taux, tauy
  real*8                            :: q0_act, emp_act, taux_act, tauy_act
  real*8                            :: tAlpha = 2e-4
  real*8                            :: sBeta = 0.*7.4e-4 
  real*8                            :: tact

  ! FIXME: Do we need these
  logical :: &
    do_mytke=.false., &
    do_mypp=.false.

  integer               :: fid    = 25
  character(len=128) :: path_data="./"

  contains
end module onedmix_variables
