subroutine drift(a,r,b)

  implicit none

  real(8), intent(in)  :: a, r(3)
  real(8), intent(out) :: b(3)

  real(8) :: ar_inv

  ar_inv = -a/dsqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
  b(:)   = r(:)*ar_inv

end subroutine drift



subroutine ave_error(x,n,ave,err)

  implicit none

  integer, intent(in)  :: n 
  real(8), intent(in)  :: x(n) 
  real(8), intent(out) :: ave, err

  real(8) :: var

  if (n.eq.1) then
    ave = x(1)
    err = 0.d0
  else
    ave = 0.d0
    ave = sum(x)/dble(n)
    var = sum((x-ave)*(x-ave))/dble(n-1)
    err = sqrt(var/dble(n))
  end if

end subroutine ave_error



subroutine random_gauss(z,n)

  implicit none

  integer, intent(in)  :: n
  real(8), intent(out) :: z(n)
  real(8)              :: u(n+1)
  real(8), parameter   :: two_pi = 2.d0*dacos(-1.d0)
  integer              :: i

  call random_number(u)

  if (iand(n,1).eq.0) then
    ! n is even
    do i = 1,n,2
      z(i)   = dsqrt(-2.d0*dlog(u(i))) 
      z(i+1) = z(i)*dsin(two_pi*u(i+1))
      z(i)   = z(i)*dcos(two_pi*u(i+1))
    end do

  else
    ! n is odd
    do i = 1,n-1,2
      z(i)   = dsqrt(-2.d0*dlog(u(i))) 
      z(i+1) = z(i)*dsin(two_pi*u(i+1))
      z(i)   = z(i)*dcos(two_pi*u(i+1))
    end do

    z(n) = dsqrt(-2.d0*dlog(u(n))) 
    z(n) = z(n)*dcos(two_pi*u(n+1))

  end if

end subroutine random_gauss
