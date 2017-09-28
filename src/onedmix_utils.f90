module onedmix_utils
  use onedmix_variables
  use onedmix_eos
  implicit none
  
  contains

!-------------------------------------------------------------------------------- 
  subroutine solve_tridiag(a,b,c,d,x,n)
  !      implicit none
  !---------------------------------------------------------------------------------
  !        a - sub-diagonal (means it is the diagonal below the main diagonal)
  !        b - the main diagonal
  !        c - sup-diagonal (means it is the diagonal above the main diagonal)
  !        d - right part
  !        x - the answer
  !        n - number of equations
  !---------------------------------------------------------------------------------
          integer,intent(in) :: n
          real*8,dimension(n),intent(in) :: a,b,c,d
          real*8,dimension(n),intent(out) :: x
          real*8,dimension(n) :: cp,dp
          real*8 :: m,fxa
          integer i
  
  ! initialize c-prime and d-prime
          cp(1) = c(1)/b(1)
          dp(1) = d(1)/b(1)
  ! solve for vectors c-prime and d-prime
           do i = 2,n
             m = b(i)-cp(i-1)*a(i)
             fxa = 1D0/m
             cp(i) = c(i)*fxa
             dp(i) = (d(i)-dp(i-1)*a(i))*fxa
           enddo
  ! initialize x
           x(n) = dp(n)
  ! solve for x from the vectors c-prime and d-prime
          do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
          end do
  end subroutine solve_tridiag

!-------------------------------------------------------------------------------- 
  subroutine interp_forcing(forc, forc_interp)
    real*8, intent(in), dimension(nforc) :: forc
    real*8, intent(out) :: forc_interp
    integer :: nf
  
    ! linear interpolation
    nf = floor(tact/force_freq)+1
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
      call calc_dens(temp(k-1), salt(k-1), pint(k), dens_km1) 
      call calc_dens(temp(k),   salt(k),   pint(k), dens_k) 
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

end module onedmix_utils
