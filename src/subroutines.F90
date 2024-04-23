!--------- Drift Vector ------------------------!

subroutine drift(a,r,nn,Rn,d)

  implicit none

  integer, intent(in)  :: nn
  real(8), intent(in)  :: a, r(3), Rn(nn,3)
  real(8), intent(out) :: d(3)

  integer :: i
  real(8) :: dist, dr(3)

  d(:) = 0.d0

  do i = 1,nn
    dr(:) = r(:) - Rn(i,:)
    dist = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

    if (dist.le.1.d-6) then
      write(*,*) '*** WARNING: drift vector diverges'
      d(:) = huge(1.d0)
    else
      d(:)   = d(:) - a*dr(:)/dist
    end if
  end do

end subroutine drift

!--------- Average and Errorbar -----------------!

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

!--------- Box-Muller ------------------------!

subroutine random_gauss(z,n)

  implicit none

  real(8), parameter   :: two_pi = 2.d0*dacos(-1.d0)

  integer, intent(in)  :: n
  real(8), intent(out) :: z(n)

  real(8) :: u(n+1)
  integer :: i

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





!!!!!!!!!!!!!!  TEST  !!!!!!!!!!!!!!!!!!!

#ifdef TEST
program test
  call test_drift
end program test
#endif

subroutine test_drift

    implicit none

    integer :: nn
    real(8) :: r(3), Rn(2,3)
    real(8) :: a = 1.d0
    real(8) :: d(3), expected_drift(3)

    nn = size(Rn, dim=1)

!--------------------------------------!

    r(:)    = (/0.d0, 0.d0, 5.d0/)
    Rn(1,:) = (/0.d0, 0.d0, 0.d0/)
    Rn(2,:) = (/0.d0, 0.d0, 1.d0/)

    expected_drift(:) = (/0.d0, 0.d0, -2.d0/)

    call drift(a, r, nn, Rn, d)

    if (d(3).ne.expected_drift(3)) stop 'Drift failed'

!--------------------------------------!

    r(:)    = (/0.d0, 0.d0, 0.d0/)
    Rn(1,:) = (/0.d0, 0.d0, 0.d0/)
    Rn(2,:) = (/0.d0, 0.d0, 1.d0/)

    expected_drift(:) = huge(-1.d0)

    call drift(a, r, nn, Rn, d)

    if (d(3).ne.expected_drift(3)) stop 'Drift failed'


    write(*,*) 'Drift ok'

end subroutine test_drift
