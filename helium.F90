#ifdef TEST_H
program test_h
  call test_potential
end program test_h
#endif


!------- POTENTIAL ENERGY ---------------------------------!
double precision function potential(r1,r2)

  implicit none

  double precision, intent(in) :: r1(3), r2(3)

  double precision :: dist1, dist2, dist12

  dist1  = dsqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3))
  dist2  = dsqrt(r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3))

  dist12 = dsqrt((r1(1)-r2(1))*(r1(1)-r2(1)) + &
                 (r1(2)-r2(2))*(r1(2)-r2(2)) + &
                 (r1(3)-r2(3))*(r1(3)-r2(3)))

  if (dist1.le.1.d-6.or.dist2.le.1.d-6) then

    potential = -huge(1.d0)
    print*, "****** Warning, potential diverges **********"

  elseif (dist12.le.1.d-6) then

    potential = huge(1.d0)
    print*, "****** Warning, potential diverges **********"

  else

    potential = 1.d0/dist12 - 2.d0/dist1 - 2.d0/dist2

  end if

end function potential


!------- KINETIC ENERGY -----------------------------------!
double precision function kinetic(a,r1,r2)

  implicit none

  double precision, intent(in) :: a, r1(3), r2(3)

  double precision :: dist1, dist2

  dist1 = dsqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3)) 
  dist2 = dsqrt(r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3)) 

  if (dist1.le.1.d-6.or.dist2.le.1.d-6) then

    write(*,*) "****** Warning, potential diverges **********"
    kinetic = huge(1.d0)

  else

    kinetic = -a*a + a/dist1 + a/dist2

  end if

end function kinetic


!------- LOCAL ENERGY -------------------------------------!
double precision function e_loc(a,r1,r2)

  implicit none

  double precision, intent(in) :: a, r1(3), r2(3)

  double precision, external :: kinetic
  double precision, external :: potential

  e_loc = kinetic(a,r1,r2) + potential(r1,r2)

end function e_loc


!------- TRIAL WAVEFUNCTION -------------------------------!
double precision function psi(a,r1,r2)

  implicit none

  double precision, intent(in) :: a, r1(3), r2(3)

  double precision :: dist1, dist2

  dist1 = dsqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3))
  dist2 = dsqrt(r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3))

  psi = dexp(-a*(dist1 + dist2))

end function psi



!subroutine test_potential
!
!    implicit none
!
!    double precision :: r1(3), r2(3)
!    double precision :: expected_output
!    double precision, external :: potential
!
!    expected_output = -1.d0/15.d0
!
!    r1(:) = (/ 2.d0, 5.d0, 14.d0 /)
!    r2(:) = (/ 1.d0, 7.d0, 4.d0 /)
!    if (potential(r) /= expected_output) stop 'Failed'
!
!    r1(:) = (/ 5.d0, 14.d0, 2.d0 /)
!    if (potential(r) /= expected_output) stop 'Failed'
!
!    r1(:) = (/ -2.d0, 5.d0, -14.d0 /)
!    if (potential(r) /= expected_output) stop 'Failed'
!
!    r1(:) = (/ 5.d0, -14.d0, -2.d0 /)
!    if (potential(r) /= expected_output) stop 'Failed'
!
!    r1(:) = (/ 0.d0, 9.d0, 12.d0 /)
!    if (potential(r) /= expected_output) stop 'Failed'
!
!    r1(:) = (/ 9.d0, -12.d0, 0.d0 /)
!    if (potential(r) /= expected_output) stop 'Failed'
!
!    r1(:) = 0.d0
!    expected_output = -huge(1.d0)
!    if (potential(r) /= expected_output) stop 'Failed r=0'
!    print *, 'potential ok'
!
!end subroutine test_potential
