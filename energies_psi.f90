!------- POTENTIAL ENERGY ---------------------------------!

real(8) function potential_1(r,Z)

  implicit none

  real(8), intent(in) :: r(3), Z

  real(8) :: dist

  dist  = dsqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))

  if (dist.le.1.d-6) then

    potential_1 = -huge(1.d0)
    write(*,*) "*** WARNING, potential diverges ***"

  else

    potential_1 = - Z/dist

  end if

end function potential_1

real(8) function potential_2(r1,r2,Z)

  implicit none

  real(8), intent(in) :: r1(3), r2(3), Z

  real(8) :: dist1, dist2, dist12

  dist1  = dsqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3))
  dist2  = dsqrt(r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3))

  dist12 = dsqrt((r1(1)-r2(1))*(r1(1)-r2(1)) + &
                 (r1(2)-r2(2))*(r1(2)-r2(2)) + &
                 (r1(3)-r2(3))*(r1(3)-r2(3)))

  if (dist1.le.1.d-6.or.dist2.le.1.d-6) then

    potential_2 = -huge(1.d0)
    write(*,*) "*** WARNING, potential diverges ***"

  elseif (dist12.le.1.d-6) then

    potential_2 = huge(1.d0)
    write(*,*) "*** WARNING, potential diverges ***"

  else

    potential_2 = 1.d0/dist12 - Z/dist1 - Z/dist2

  end if

end function potential_2


!------- KINETIC ENERGY -----------------------------------!

real(8) function kinetic_1(a,r)

  implicit none

  real(8), intent(in) :: a, r(3)

  real(8) :: dist

  dist = dsqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)) 

  if (dist.le.1.d-6) then

    write(*,*) "*** WARNING, potential diverges ***"
    kinetic_1 = huge(1.d0)

  else

    kinetic_1 = a/dist - a*a*0.5d0

  end if

end function kinetic_1

real(8) function kinetic_2(a,r1,r2)

  implicit none

  real(8), intent(in) :: a, r1(3), r2(3)

  real(8) :: dist1, dist2

  dist1 = dsqrt(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3)) 
  dist2 = dsqrt(r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3)) 

  if (dist1.le.1.d-6.or.dist2.le.1.d-6) then

    write(*,*) "*** WARNING, potential diverges ***"
    kinetic_2 = huge(1.d0)

  else

    kinetic_2 = -a*a + a/dist1 + a/dist2

  end if

end function kinetic_2


!------- LOCAL ENERGY -------------------------------------!

real(8) function e_loc_1(a,r,Z)

  implicit none

  real(8), intent(in) :: a, r(3), Z

  real(8), external :: kinetic_1
  real(8), external :: potential_1

  e_loc_1 = kinetic_1(a,r) + potential_1(r,Z)

end function e_loc_1

real(8) function e_loc_2(a,r1,r2,Z)

  implicit none

  real(8), intent(in) :: a, r1(3), r2(3), Z

  real(8), external :: kinetic_2
  real(8), external :: potential_2

  e_loc_2 = kinetic_2(a,r1,r2) + potential_2(r1,r2,Z)

end function e_loc_2


!------- TRIAL WAVEFUNCTION -------------------------------!

real(8) function psi_1(a,r,nn,Rn)

  implicit none

  integer, intent(in) :: nn
  real(8), intent(in) :: a, r(3), Rn(nn,3)

  integer :: i
  real(8) :: dist, dr(3)

  psi_1 = 0.d0

  do i = 1,nn
    dr(:) = r(:) - Rn(i,:)

    dist = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
    psi_1 = psi_1 + dexp(-a*dist)
  end do

end function psi_1

real(8) function psi_2(a,r1,r2,nn,Rn)

  implicit none

  integer, intent(in) :: nn
  real(8), intent(in) :: a, r1(3), r2(3), Rn(nn,3)

  integer :: i
  real(8) :: dist1, dist2, dr1(3), dr2(3)

  psi_2 = 0.d0

  do i = 1,nn
    dr1(:) = r1(:) - Rn(i,:) 
    dr2(:) = r2(:) - Rn(i,:)

    dist1 = dsqrt(dr1(1)*dr1(1) + dr1(2)*dr1(2) + dr1(3)*dr1(3))
    dist2 = dsqrt(dr2(1)*dr2(1) + dr2(2)*dr2(2) + dr2(3)*dr2(3))

    psi_2 = psi_2 + dexp(-a*(dist1 + dist2))
  end do

end function psi_2
