program x2t2025
  ! x is a/a0 = a version 2025
  implicit none
  integer(kind=4)::i
  integer(kind=4),parameter::na=10**8
  real(kind=8)::a,da,integral
  real(kind=8),dimension(1:na)::t
  real(kind=8),parameter::TH=14.429e9 ! TH = 1/H_0 (y)
  real(kind=8),parameter::OmL=0.6889, Omm=0.3071, Omr=9.0824e-5
  real(kind=8),parameter::OmK=0.0000, half = 0.5d0, three = 3.0d0
  real(kind=8),parameter::zero=0.0d0, one=1.0d0, tenth=0.1d0
  real(kind=8),parameter::w0 = - 0.752d0, wa = - 0.86d0 ! DESI
  integral = zero
  da = one/na; a = half*da
  do i = 1, na
     integral = integral + da*TH/sqrt(OmL*(a**(-three*(one + w0 + wa))*exp(-three*wa*(one - a))*a**2 + OmK + Omm/a + Omr/a**2))
     t(i) = integral ! time in years to get to a/a_0 = x
     a = a + da
  end do
   open(7,file='desi.dat')
   do i = 1, na, 1000
  !                    a/a_0     t      z(t)                              T_r(t)
     write(7,*) (half + i - 1)*da, t(i), one/((i-half)*da) - one, 2.7255/((i+1)*da)
  end do
  close(7)
end program x2t2025
