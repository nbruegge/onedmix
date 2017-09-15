module onedmix_cvmix_tke
  use cvmix_tke
  use onedmix_variables
  use onedmix_eos
  implicit none

REAL*8 :: forc_tke_surf, bottom_fric, forc_rho_surf
REAL*8, DIMENSION(:), ALLOCATABLE, TARGET :: &
  tke_diss                                        ,&
  tke_iw_forcing                                  ,&
  tke_iwe                                         ,&
  tke_iw_alpha_c                                  ,&
  cvmix_dummy_1                                   ,&
  cvmix_dummy_2                                   ,&
  cvmix_dummy_3                                   ,&
  ! nils tke diagnostics
  tke                                             ,&
  tke_Lmix                                        ,&
  tke_Pr                                          ,&
  tke_Tbpr                                        ,&
  tke_Tspr                                        ,&
  tke_Tdif                                        ,&
  tke_Tdis                                        ,&
  tke_Twin                                        ,&
  tke_Tiwf                                        ,&
  tke_Tbck                                        ,&
  tke_Ttot                                        ,&
  tke_avo                                         ,&
  tke_dvo                                         
  ! end nils tke diagnostics

! VMIX_NL/CVMIX_TKE_PARAM namelist parameters
REAL*8 :: &
  c_k, c_eps, alpha_tke, mxl_min, kappaM_min, kappaM_max, cd, tke_surf_min, tke_min
INTEGER  :: &
  tke_mxl_choice
LOGICAL :: &
  l_tke_active, only_tke, use_ubound_dirichlet, use_lbound_dirichlet


  contains
!-------------------------------------------------------------------------------- 
  subroutine setup_cvmix_tke
    character(len=128)    :: fname
    namelist /tke_paras/                                                   &
        c_k, c_eps, alpha_tke, mxl_min, kappaM_min, kappaM_max, tke_mxl_choice   &
      , cd, tke_surf_min, tke_min, l_tke_active                                  &
      , only_tke, use_ubound_dirichlet, use_lbound_dirichlet
    open(fid, file="./onedmix.nl", status="old", action='read')
    read(fid, nml=tke_paras)
    close(fid)

    ! if IDEMIX is active these variables are calculated in cvmix_idemix 
    ! otherwise they stay zero
    ALLOCATE(tke_iw_forcing(nz+1))
      tke_iw_forcing(:)=0.0
    ALLOCATE(tke_iwe(nz+1))
      tke_iwe(:)=0.0
    ALLOCATE(tke_iw_alpha_c(nz+1))
      tke_iw_alpha_c(:)=0.0
   
    ALLOCATE(tke_diss(nz+1))
      tke_diss(:)=0.0
    ALLOCATE(cvmix_dummy_1(nz+1))
     cvmix_dummy_1(:)=0.0
    ALLOCATE(cvmix_dummy_2(nz+1))
     cvmix_dummy_2(:)=0.0
    ALLOCATE(cvmix_dummy_3(nz+1))
     cvmix_dummy_3(:)=0.0
   
    ALLOCATE(tke(nz+1))
     tke(:)=0.0
    ALLOCATE(tke_Lmix(nz+1))
     tke_Lmix(:)=0.0
    ALLOCATE(tke_Pr(nz+1))
     tke_Pr(:)=0.0
    ALLOCATE(tke_Tbpr(nz+1))
     tke_Tbpr(:)=0.0
    ALLOCATE(tke_Tspr(nz+1))
     tke_Tspr(:)=0.0
    ALLOCATE(tke_Tdif(nz+1))
     tke_Tdif(:)=0.0
    ALLOCATE(tke_Tdis(nz+1))
     tke_Tdis(:)=0.0
    ALLOCATE(tke_Twin(nz+1))
     tke_Twin(:)=0.0
    ALLOCATE(tke_Tiwf(nz+1))
     tke_Tiwf(:)=0.0
    ALLOCATE(tke_Tbck(nz+1))
     tke_Tbck(:)=0.0
    ALLOCATE(tke_Ttot(nz+1))
     tke_Ttot(:)=0.0
    ALLOCATE(tke_avo(nz+1))
     tke_avo(:)=0.0
    ALLOCATE(tke_dvo(nz+1))
     tke_dvo(:)=0.0
   
    !tke forcing fields
    forc_tke_surf = 0.0
    forc_rho_surf = 0.0
    bottom_fric   = 0.0

    ! read initial tke
    fname = trim(path_data) // "etke0.txt"
    open(fid, file=fname, status="old", action='read')
    read(fid, *) tke
    close(fid)
    
    CALL init_tke(c_k            = c_k,            &
                  c_eps          = c_eps,          &
                  cd             = cd,             &
                  alpha_tke      = alpha_tke,      &
                  mxl_min        = mxl_min,        &
                  kappaM_min     = kappaM_min,     &
                  kappaM_max     = kappaM_max,     &
                  tke_mxl_choice = tke_mxl_choice, &
                  use_ubound_dirichlet = use_ubound_dirichlet, &
                  use_lbound_dirichlet = use_lbound_dirichlet, &
                  only_tke       = only_tke,       &
                  tke_min        = tke_min,        &
                  tke_surf_min   = tke_surf_min    )
  end subroutine setup_cvmix_tke

!-------------------------------------------------------------------------------- 
  subroutine calc_cvmix_tke
    REAL*8, DIMENSION(nz+1) :: &
      avo_old, &
      dvo_old!, &
    integer :: i, j, tstep_count
    i = 1
    j = 1
    tstep_count = 1
    avo_old = 0.0
    dvo_old = 0.0
    
    !tke forcing fields
    forc_tke_surf = cd * sqrt( (taux_act**2 + tauy_act**2)**(3./2.)  )
    forc_rho_surf = 0.0
    bottom_fric   = 0.0

    ! main cvmix call to calculate tke
    CALL cvmix_coeffs_tke( &
                           i = i,                          &
                           j = j,                          &
                           tstep_count = tstep_count,      &
                           tke_diss_out = tke_diss,        & ! (inout)
                           tke_out      = tke,             & ! (inout)
                           ! FIXME: nils: exchange cvmix_dummy_1 with avo later
                           KappaM_out   = tke_avo,         & ! (inout)
                           KappaH_out   = tke_dvo,         & ! (inout)
                           !KappaM_out   = avo,             & ! (inout)
                           !KappaH_out   = dvo,             & ! (inout)
                           cvmix_int_1  = cvmix_dummy_1,   & ! (out)
                           cvmix_int_2  = cvmix_dummy_2,   & ! (out)
                           cvmix_int_3  = cvmix_dummy_3,   & ! (out)
                           dzw          = dzw,             &
                           dzt          = dzt,             &
                           nlev         = nz,              &
                           max_nlev     = nz,              &
                           old_tke      = tke,             &
                           old_tke_diss = tke_diss,        &
                           Ssqr         = S2,              &
                           Nsqr         = N2,              &
                           tke_Tbpr     = tke_Tbpr,        &
                           tke_Tspr     = tke_Tspr,        &
                           tke_Tdif     = tke_Tdif,        &
                           tke_Tdis     = tke_Tdis,        &
                           tke_Twin     = tke_Twin,        &
                           tke_Tiwf     = tke_Tiwf,        &
                           tke_Tbck     = tke_Tbck,        &
                           tke_Ttot     = tke_Ttot,        &
                           tke          = tke,             &
                           tke_Lmix     = tke_Lmix,        &
                           tke_Pr       = tke_Pr,          &
                           forc_tke_surf= forc_tke_surf,   &
                           ! FIXME: nils: needs to be set some where
                           E_iw         = tke_iwe,         &
                           dtime        = dt,              &
                           bottom_fric  = bottom_fric,     &
                           old_kappaM   = avo_old,         &
                           old_KappaH   = dvo_old,         & 
                           iw_diss      = tke_iw_forcing,  & 
                           forc_rho_surf= forc_rho_surf,   &
                           !Kappa_GM     = Kappa_GM,             & ! FIXME: optional
                           rho_ref      = rho0,            &
                           grav         = grav,            &
                           alpha_c      = tke_iw_alpha_c   &
                         )
  end subroutine calc_cvmix_tke

!-------------------------------------------------------------------------------- 
  subroutine write_snap_cvmix_tke
  end subroutine write_snap_cvmix_tke

end module onedmix_cvmix_tke
