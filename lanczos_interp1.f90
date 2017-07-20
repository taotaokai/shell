subroutine lanczos_interp1(x, n, dt, ti, ni, na, xi)
!   Ideal bandlimited (Sinc) interpolation.
!   
!   :type dt: float
!   :param dt: sampling interval.
!   :type x: float
!   :param x: 1d array of samples
!   :type ti: float
!   :param ti: 1d array of interpolation times. 
!       The zero time is at the first sample of x.
!   :type na: int
!   :param na: half width of sinc window (number of samples)
!       the larger this value the better reconstruction, but the longer
!       duration(~na*dt) of the artefacts at the beginning and end of signal.
!   :rtype: float
!   :returns: 1d array of interpolation results

implicit none

! arguments
integer :: n, ni, na
!f2py integer intent(hide),depend(x) :: n = len(x)
!f2py integer intent(hide),depend(ti) :: ni = len(ti)
!f2py integer optional,intent(in) :: na = 10
real(kind=8), intent(in) :: x(0:n-1)
real(kind=8), intent(in) :: dt
real(kind=8), intent(in) :: ti(0:ni-1)
real(kind=8), intent(out) :: xi(0:ni-1)

integer, parameter :: dp = kind(0.d0)
real(dp), parameter :: PI = 4 * atan(1.0_dp)
integer :: i, j, i0, i1
real(dp) :: ii, s
real(dp) :: PI_s, PI_s_na, sinc_s, sinc_s_na

xi = 0.0_dp
do i = 0, ni-1
    ii = ti(i)/dt
    ! index range of sinc window
    i0 = max(floor(ii)-na+1, 0)
    i1 = min(floor(ii)+na, n-1)
    do j = i0, i1
        ! sinc(s): bandlimited pulse response at ti(i) from x(j)
        s = ii - j
        PI_s = PI * s
        sinc_s = 1.0_dp
        if (PI_s .ne. 0.0_dp) sinc_s = sin(PI_s)/PI_s
        ! sinc(s/na): window around ti(i)
        PI_s_na = PI_s / na
        sinc_s_na = 1.0_dp
        if (PI_s_na .ne. 0.0_dp) sinc_s_na = sin(PI_s_na)/PI_s_na
        ! reconstructed signal at ti(i)
        xi(i) = xi(i) + sinc_s*sinc_s_na*x(j)
    enddo
enddo

end subroutine