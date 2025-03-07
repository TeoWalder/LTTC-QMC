!------- POTENTIAL ENERGY ---------------------------------!

real(8) function potential(r,ne,Rn,nn,Z)

  implicit none

  integer, intent(in) :: ne, nn
  real(8), intent(in) :: r(ne,3), Rn(nn,3), Z(nn)

  integer :: i, j
  real(8) :: dist, dr(3)

  potential = 0.d0

  ! Electron-Nucleus interaction
  do i = 1,nn
    do j = 1,ne
      dr(:) = r(j,:) - Rn(i,:)
      dist  = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

      if (dist.le.1.d-6) then
        potential = potential - huge(1.d0)
        write(*,*) '*** WARNING: el-Nu potential diverges'
      else
        potential = potential - Z(i)/dist
      end if

    end do
  end do

  ! Electron-Electron interaction
  do i = 1,ne-1
    do j = i+1,ne
      dr(:) = r(j,:) - r(i,:)
      dist  = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

      if (dist.le.1.d-6) then
        potential = potential + huge(1.d0)
        write(*,*) '*** WARNING: el-el potential diverges'
      else
        potential = potential + 1.d0/dist
      end if

    end do
  end do

  ! Nucleus-Nucleus interaction
  do i = 1,nn-1
    do j = i+1,nn
      dr(:) = Rn(j,:) - Rn(i,:)
      dist  = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

      if (dist.le.1.d-6) then
        potential = potential + huge(1.d0)
        write(*,*) '*** WARNING: Nu-Nu potential diverges'
      else
        potential = potential + Z(i)*Z(i)/dist
      end if

    end do
  end do

end function potential


!------- KINETIC ENERGY -----------------------------------!

real(8) function kinetic(a,r,ne,Rn,nn)

  implicit none

  integer, intent(in) :: ne, nn
  real(8), intent(in) :: a, r(ne,3), Rn(nn,3)

  integer :: i, j
  real(8) :: dist, dr(3), k1, k2
  real(8) :: psi_tot, psi_Ri, ratiopsi

  real(8), external :: psi, psi_Rn

  ! constant term
  k1 = -a*a*0.5d0 * dble(ne)
  
  ! psi dependent term
  k2 = 0.d0
  psi_tot = psi(a,r,nn,ne,Rn)

  if (psi_tot.le.1.d-6) then
    write(*,*) '*** WARNING: kinetic energy diverges'
    k2 = huge(1.d0)
  end if

  do i = 1,nn
    psi_Ri = psi_Rn(a,r,ne,Rn(i,:))
    ratiopsi = psi_Ri/psi_tot

    do j = 1,ne
      dr(:) = r(j,:) - Rn(i,:)
      dist = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
  
      if (dist.le.1.d-6) then
        write(*,*) '*** WARNING: kinetic energy diverges'
        k2 = k2 + huge(1.d0)
      else
        k2 = k2 + a/dist * ratiopsi
      end if
    end do
    
  end do

  kinetic = k1 + k2

end function kinetic


!------- LOCAL ENERGY -------------------------------------!

real(8) function e_loc(a,r,ne,Rn,nn,Z)

  implicit none

  integer, intent(in) :: ne, nn
  real(8), intent(in) :: a, r(ne,3), Rn(nn,3), Z(nn)

  real(8), external :: kinetic
  real(8), external :: potential

  e_loc = kinetic(a,r,ne,Rn,nn) + potential(r,ne,Rn,nn,Z)

end function e_loc


!------- IONIC WAVEFUNCTION -------------------------------!

real(8) function psi(a,r,nn,ne,Rn)

  implicit none

  integer, intent(in) :: nn, ne
  real(8), intent(in) :: a, r(ne,3), Rn(nn,3)

  integer :: i
  real(8) :: hp

  real(8), external :: psi_Rn

  psi = 0.d0

  do i = 1,nn
    hp  = psi_Rn(a,r,ne,Rn(i,:))
    psi = psi + hp
  end do

end function psi

! hartree product per nucleus
real(8) function psi_Rn(a,r,ne,Rn)

  implicit none

  integer, intent(in) :: ne
  real(8), intent(in) :: a, r(ne,3), Rn(3)

  integer :: i
  real(8) :: dist, disti, dr(3)

  dist = 0.d0

  do i = 1,ne
    dr(:) = r(i,:) - Rn(:)
    disti = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
    dist = dist + disti
  end do

  psi_Rn = dexp(-a*dist)

end function psi_Rn



!!!!!!!!!!!!!!  TEST  !!!!!!!!!!!!!!!!!!!

#ifdef TEST
program test
  call test_energies
end program test
#endif

subroutine test_energies

    implicit none

    integer :: ne
    integer :: nn
    real(8) :: a = 1.d0
    real(8) :: r(2,3), Rn(2,3), Z(2)
    real(8) :: expected_pot, expected_kin

    real(8), external :: potential, kinetic, psi, psi_Rn

    ne = size(r, dim=1)
    nn = size(Rn, dim=1)

    Z(:)    = (/1.d0, 1.d0/)
    r(1,:)  = (/0.d0, 0.d0, 2.d0/)
    r(2,:)  = (/0.d0, 0.d0, -1.d0/)
    Rn(1,:) = (/0.d0, 0.d0, 0.d0/)
    Rn(2,:) = (/0.d0, 0.d0, 1.d0/)

    expected_pot = -1.6666666666666665d0
    expected_kin = 0.5d0

    if (potential(r,ne,Rn,nn,Z).ne.expected_pot) stop 'Potential Failed'
    if (kinetic(a,r,ne,Rn,nn).ne.expected_kin)   stop 'Kinetic Failed'

!------------------------------------------------!

    Z(:)    = (/1.d0, 1.d0/)
    r(1,:)  = (/0.d0, 0.d0, -1.d0/)
    r(2,:)  = (/0.d0, 0.d0, -1.d0/)
    Rn(1,:) = (/0.d0, 0.d0, 0.d0/)
    Rn(2,:) = (/0.d0, 0.d0, 1.d0/)

    expected_pot = huge(1.d0)
    expected_kin = 0.88079707797788243d0

    if (potential(r,ne,Rn,nn,Z).ne.expected_pot) stop 'Potential Failed'
    if (kinetic(a,r,ne,Rn,nn).ne.expected_kin)   stop 'Kinetic Failed'

    write(*,*) 'Energies ok'

end subroutine test_energies
