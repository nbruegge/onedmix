module onedmix_eos
  implicit none

  contains

!-------------------------------------------------------------------------------- 
  subroutine calc_dens(itemp, isalt, pres0, odens, inz)
    !use onedmix_setup
    !use eos
    integer, intent(in)                 :: inz 
    real*8, intent(in), dimension(inz)  :: itemp, isalt 
    real*8, intent(in)                  :: pres0
    real*8, intent(out), dimension(inz) :: odens
    integer :: k
    ! FIXME: linear eos for the moment
    !do k=1,inz
    !  odens(k) = rho0 * ( 1.d0 - tAlpha*itemp(k) + sBeta*isalt(k) )
    !enddo

    do k=1,inz
      call potrho(itemp(k), isalt(k), pres0, odens(k))
    end do

    !pbcl(1) = 0.0
    !do k=2,inz
    !  pbcl(k) = pbcl(k-1) + grav*0.5*(odens(k-1)+odens(k))*dzt(k)
    !end do
  end subroutine calc_dens

!-------------------------------------------------------------------------------- 
! from MPIOM
! compute density from potential temperature directly
  SUBROUTINE potrho(tpot, sal, p, rho)
    !INTEGER, INTENT(in) :: nz 
    REAL*8, INTENT(in) :: tpot, sal, p
    REAL*8, INTENT(out) :: rho
    REAL*8, PARAMETER :: &
         a_a1=3.6504E-4, a_a2=8.3198E-5, a_a3=5.4065E-7, &
         a_a4=4.0274E-9, &
         a_b1=1.7439E-5, a_b2=2.9778E-7, &
         a_c1=8.9309E-7, a_c2=3.1628E-8, a_c3=2.1987E-10, &
         a_d=4.1057E-9, &
         a_e1=1.6056E-10, a_e2=5.0484E-12

    REAL*8, PARAMETER :: &
         r_a0=999.842594, r_a1=6.793952e-2, r_a2=-9.095290e-3, &
         r_a3=1.001685e-4, r_a4=-1.120083e-6, r_a5=6.536332e-9, &
         r_b0=8.24493e-1, r_b1=-4.0899e-3, r_b2=7.6438e-5, &
         r_b3=-8.2467e-7, r_b4=5.3875e-9, &
         r_c0=-5.72466e-3, r_c1=1.0227e-4, r_c2=-1.6546e-6, &
         r_d0=4.8314e-4, &
         r_e0=19652.21, r_e1=148.4206, r_e2=-2.327105, &
         r_e3=1.360477e-2, r_e4=-5.155288e-5, &
         r_f0=54.6746, r_f1=-0.603459, r_f2=1.09987e-2, &
         r_f3=-6.1670e-5, &
         r_g0=7.944e-2, r_g1=1.6483e-2, r_g2=-5.3009e-4, &
         r_h0=3.239908, r_h1=1.43713e-3, r_h2=1.16092e-4, &
         r_h3=-5.77905e-7, &
         r_ai0=2.2838e-3, r_ai1=-1.0981e-5, r_ai2=-1.6078e-6, &
         r_aj0=1.91075e-4, &
         r_ak0=8.50935e-5, r_ak1=-6.12293e-6, r_ak2=5.2787e-8, &
         r_am0=-9.9348e-7, r_am1=2.0816e-8, r_am2=9.1697e-10

    REAL*8 :: dc, dv, dvs, fne, fst, qc, qn3, qnq, qv, qvs, &
         s, s3h, t, tpo
    INTEGER :: k

    qc = p * (a_a1 + p * (a_c1 - a_e1 * p))
    qv = p * (a_b1 - a_d * p)
    dc = 1. + p * (-a_a2 + p * (a_c2 - a_e2 * p))
    dv = a_b2 * p
    qnq  = -p * (-a_a3 + p * a_c3)
    qn3  = -p * a_a4

    tpo = tpot
    qvs = qv*(sal - 35.0) + qc
    dvs = dv*(sal - 35.0) + dc
    t   = (tpo + qvs)/dvs
    fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpo
    fst = dvs + t*(2.0*qnq + 3.*qn3*t)
    t = t - fne/fst
    s = MAX(sal, 0.0)
    s3h=SQRT(s**3)

    rho = &
         (r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))&
         & + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))    &
         & + r_d0 * s**2                                                     &
         & + s3h * (r_c0 + t * (r_c1 + r_c2 * t)))                           &
         / (1.                                                            &
         &  - p / (p * (r_h0 + t * (r_h1 + t * (r_h2 + t * r_h3))            &
         &              + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))              &
         &              + r_aj0 * s3h                                        &
         &              + (r_ak0 + t * (r_ak1 + t * r_ak2)                   &
         &              + s * (r_am0 + t * (r_am1 + t * r_am2))) * p)        &
         &         + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))  &
         &         + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))         &
         &         + s3h * (r_g0 + t * (r_g1 + r_g2 * t))))
  END SUBROUTINE potrho
end module onedmix_eos
