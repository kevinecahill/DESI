program x2t2025
  ! x is a/a0 = a version 2025
  implicit none
  integer(kind=4)::i
  integer(kind=4),parameter::nx=10**8
  real(kind=8)::x,dx,integral
  real(kind=8),dimension(1:nx)::t
  real(kind=8),parameter::TH=14.429e9 ! TH = 1/H_0 (y)
  real(kind=8),parameter::OmL=0.6889, Omm=0.3071, Omr=9.0824e-5
  real(kind=8),parameter::OmK=0.0000, half = 0.5d0
  real(kind=8),parameter::zero=0.0d0, one=1.0d0, tenth=0.1d0
  integral = zero
  dx = one/nx; x = half*dx
  do i = 1, nx
     integral = integral + dx*TH/sqrt(OmL*x**2 + OmK + Omm/x + Omr/x**2)
     t(i) = integral ! time in years to get to a/a_0 = x
     x = x + dx
  end do
   open(7,file='Planck.dat')
   do i = 1, nx, 1000
  !                    a/a_0     t      z(t)                              T_r(t)
     write(7,*) (half + i - 1)*dx, t(i), one/((i-half)*dx) - one, 2.7255/((i+1)*dx)
  end do
  close(7)
end program x2t2025
