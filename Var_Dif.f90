!----- Variational Monte Carlo ----------------------------!

subroutine VMC(a,nmax,dt,energy,accep,Rn,ne,nn,Z)

  implicit none

  integer   , intent(in)  :: ne, nn
  real(8)   , intent(in)  :: a, dt, Rn(nn,3), Z(nn)
  integer(8), intent(in)  :: nmax 
  real(8)   , intent(out) :: energy, accep

  integer(8) :: istep
  integer(8) :: n_accep
  real(8)    :: sq_dt, chi(ne,3), dq_old(ne), prod(ne), u
  real(8)    :: psi_old, psi_new, dq_new(ne), argexpo, q
  real(8)    :: r_old(ne,3), r_new(ne,3)
  real(8)    :: d_old(ne,3), d_new(ne,3)

  real(8), external :: e_loc, psi

  integer :: i

  sq_dt = dsqrt(dt)

  ! Initialization

  energy = 0.d0
  n_accep = 0_8

  do i = 1,ne
    call random_gauss(r_old(i,:), 3)

    call drift(a, r_old(i,:), nn, Rn, d_old(i,:))

    dq_old(i) = d_old(i,1)*d_old(i,1) + &
                d_old(i,2)*d_old(i,2) + &
                d_old(i,3)*d_old(i,3)
  end do

  psi_old = psi(a, r_old, nn, ne, Rn, Z)

  ! Propagation
  do istep = 1,nmax

    energy = energy + e_loc(a,r_old,ne,Rn,nn,Z)

    do i = 1,ne
      call random_gauss(chi(i,:), 3)
      r_new(i,:) = r_old(i,:) + dt*d_old(i,:) + chi(i,:)*sq_dt
  
      call drift(a, r_new(i,:), nn, Rn, d_new(i,:))

      dq_new(i) = d_new(i,1)*d_new(i,1) + &
                  d_new(i,2)*d_new(i,2) + &
                  d_new(i,3)*d_new(i,3)
    end do

    psi_new = psi(a, r_new, nn, ne, Rn, Z)

    ! Metropolis
    prod(:) = (d_new(:,1) + d_old(:,1))*(r_new(:,1) - r_old(:,1)) + &
              (d_new(:,2) + d_old(:,2))*(r_new(:,2) - r_old(:,2)) + &
              (d_new(:,3) + d_old(:,3))*(r_new(:,3) - r_old(:,3))

    argexpo = 0.5d0*sum((dq_new(:) - dq_old(:))*dt + prod(:))

    q = psi_new/psi_old
    q = dexp(-argexpo)*q*q

    call random_number(u)

    if (u.le.q) then
      n_accep = n_accep + 1_8

      r_old(:,:) = r_new(:,:)
      d_old(:,:) = d_new(:,:)
      dq_old(:)  = dq_new(:)
      psi_old    = psi_new

    end if

  end do

  energy = energy/dble(nmax)
  accep  = dble(n_accep)/dble(nmax)

end subroutine VMC


!----- Pure Diffusion Monte Carlo -------------------------!

subroutine PDMC(a,dt,nmax,energy,accep,tau,E_ref,Rn,ne,nn,Z)

  implicit none

  integer   , intent(in)  :: ne, nn
  real(8)   , intent(in)  :: a, dt, tau, E_ref, Rn(nn,3), Z
  integer(8), intent(in)  :: nmax 
  real(8)   , intent(out) :: energy, accep

  integer(8) :: istep
  integer(8) :: n_accep
  real(8)    :: sq_dt, chi(ne,3), dq_old(ne), prod(ne), u
  real(8)    :: psi_old, psi_new, dq_new(ne), argexpo, q
  real(8)    :: r_old(ne,3), r_new(ne,3)
  real(8)    :: d_old(ne,3), d_new(ne,3)
  real(8)    :: e, w, normalization, tau_current

  real(8), external :: e_loc, psi

  integer :: i

  sq_dt = dsqrt(dt)

  ! Initialization
  energy  = 0.d0
  n_accep = 0_8

  normalization = 0.d0
  w             = 1.d0
  tau_current   = 0.d0

  do i = 1,ne
    call random_gauss(r_old(i,:), 3)

    call drift(a, r_old(i,:), nn, Rn, d_old(i,:))

    dq_old(i) = d_old(i,1)*d_old(i,1) + &
                d_old(i,2)*d_old(i,2) + &
                d_old(i,3)*d_old(i,3)
  end do

  psi_old = psi(a, r_old, nn, ne, Rn, Z)

  ! Propagation
  do istep = 1,nmax

    e = e_loc(a,r_old,ne,Rn,nn,Z)
    w = w*dexp(-dt*(e - E_ref))

    normalization = normalization + w
    energy = energy + w*e

    tau_current = tau_current + dt

    ! Reset when tau is reached
    if (tau_current.gt.tau) then
      w           = 1.d0
      tau_current = 0.d0
    endif

    do i = 1,ne
      call random_gauss(chi(i,:), 3)
      r_new(i,:) = r_old(i,:) + dt*d_old(i,:) + chi(i,:)*sq_dt
  
      call drift(a, r_new(i,:), nn, Rn, d_new(i,:))
      dq_new(i) = d_new(i,1)*d_new(i,1) + &
                  d_new(i,2)*d_new(i,2) + &
                  d_new(i,3)*d_new(i,3)
    end do

    psi_new = psi(a, r_new, nn, ne, Rn, Z)

    ! Metropolis
    prod(:) = (d_new(:,1) + d_old(:,1))*(r_new(:,1) - r_old(:,1)) + &
              (d_new(:,2) + d_old(:,2))*(r_new(:,2) - r_old(:,2)) + &
              (d_new(:,3) + d_old(:,3))*(r_new(:,3) - r_old(:,3))

    argexpo = 0.5d0*sum((dq_new(:) - dq_old(:))*dt + prod(:))

    q = psi_new/psi_old
    q = dexp(-argexpo)*q*q

    call random_number(u)

    if (u.le.q) then

      n_accep = n_accep + 1_8

      r_old(:,:) = r_new(:,:)
      d_old(:,:) = d_new(:,:)
      dq_old(:)  = dq_new(:)
      psi_old    = psi_new

    end if

  end do

  energy = energy/normalization
  accep  = dble(n_accep)/dble(nmax)

end subroutine PDMC
