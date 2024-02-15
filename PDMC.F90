subroutine pdmc(a, dt, nmax, energy, accep, tau, E_ref)

  implicit none

  real(8)   , intent(in)  :: a, dt, tau
  integer(8), intent(in)  :: nmax 
  real(8)   , intent(out) :: energy, accep
  real(8)   , intent(in)  :: E_ref

  integer(8) :: istep
  integer(8) :: n_accep
  real(8)    :: sq_dt, chi(3), xi(3), dq1_old, dq2_old, prod1, prod2, u
  real(8)    :: psi_old, psi_new, dq1_new, dq2_new, argexpo1, argexpo2, q
  real(8)    :: r1_old(3), r2_old(3), r1_new(3), r2_new(3)
  real(8)    :: d1_old(3), d1_new(3), d2_old(3), d2_new(3)
  real(8)    :: e, w, normalization, tau_current

  real(8), external :: e_loc, psi

  sq_dt = dsqrt(dt)

  ! Initialization
  energy  = 0.d0
  n_accep = 0_8
  normalization = 0.d0

  w           = 1.d0
  tau_current = 0.d0

!  call random_gauss(r_old,3)
!
!  call drift(a,r_old,d_old)
!  d2_old  = d_old(1)*d_old(1) + &
!            d_old(2)*d_old(2) + &
!            d_old(3)*d_old(3)
!
!  psi_old = psi(a,r_old)

  call random_gauss(r1_old, 3)
  call random_gauss(r2_old, 3)

  call drift(a, r1_old, d1_old)
  call drift(a, r2_old, d2_old)

  dq1_old = d1_old(1)*d1_old(1) + &
            d1_old(2)*d1_old(2) + &
            d1_old(3)*d1_old(3)
  dq2_old = d2_old(1)*d2_old(1) + &
            d2_old(2)*d2_old(2) + &
            d2_old(3)*d2_old(3)

  psi_old = psi(a, r1_old, r2_old)



!-------------- DIFFUSION --------------!

  do istep = 1,nmax

     e = e_loc(a, r1_old, r2_old)
     w = w * dexp(-dt*(e - E_ref))

     normalization = normalization + w
     energy = energy + w*e

     tau_current = tau_current + dt

     ! Reset when tau is reached
     if (tau_current.gt.tau) then
        w           = 1.d0
        tau_current = 0.d0
     endif

     call random_gauss(chi, 3)
     call random_gauss( xi, 3)
     r1_new(:) = r1_old(:) + dt*d1_old(:) + chi(:)*sq_dt
     r2_new(:) = r2_old(:) + dt*d2_old(:) +  xi(:)*sq_dt

     call drift(a, r1_new, d1_new)
     call drift(a, r2_new, d2_new)
     dq1_new = d1_new(1)*d1_new(1) + &
               d1_new(2)*d1_new(2) + &
               d1_new(3)*d1_new(3)
     dq2_new = d2_new(1)*d2_new(1) + &
               d2_new(2)*d2_new(2) + &
               d2_new(3)*d2_new(3)

     psi_new = psi(a, r1_new, r2_new)

!     call random_gauss(chi,3)
!     r_new(:) = r_old(:) + dt*d_old(:) + chi(:)*sq_dt
!
!     call drift(a,r_new,d_new)
!     d2_new = d_new(1)*d_new(1) + &
!              d_new(2)*d_new(2) + &
!              d_new(3)*d_new(3)
!
!     psi_new = psi(a,r_new)

     !----- Metropolis -------------------------------------!
     prod1 = (d1_new(1) + d1_old(1))*(r1_new(1) - r1_old(1)) + &
             (d1_new(2) + d1_old(2))*(r1_new(2) - r1_old(2)) + &
             (d1_new(3) + d1_old(3))*(r1_new(3) - r1_old(3))
     prod2 = (d2_new(1) + d2_old(1))*(r2_new(1) - r2_old(1)) + &
             (d2_new(2) + d2_old(2))*(r2_new(2) - r2_old(2)) + &
             (d2_new(3) + d2_old(3))*(r2_new(3) - r2_old(3))

     argexpo1 = 0.5d0*(dq1_new - dq1_old)*dt + prod1
     argexpo2 = 0.5d0*(dq2_new - dq2_old)*dt + prod2

     q = psi_new/psi_old
     q = dexp(-argexpo1 -argexpo2)*q*q

     call random_number(u)

     if (u.le.q) then

        n_accep = n_accep + 1_8

        r1_old(:) = r1_new(:)
        r2_old(:) = r2_new(:)
        d1_old(:) = d1_new(:)
        d2_old(:) = d2_new(:)
        dq1_old   = dq1_new
        psi_old   = psi_new

     end if

!     prod = (d_new(1) + d_old(1))*(r_new(1) - r_old(1)) + &
!            (d_new(2) + d_old(2))*(r_new(2) - r_old(2)) + &
!            (d_new(3) + d_old(3))*(r_new(3) - r_old(3))
!
!     argexpo = 0.5d0*(d2_new - d2_old)*dt + prod
!
!     q = psi_new/psi_old
!     q = dexp(-argexpo)*q*q
!
!     call random_number(u)
!
!     if (u.le.q) then
!
!        n_accep = n_accep + 1_8
!
!        r_old(:) = r_new(:)
!        d_old(:) = d_new(:)
!        d2_old   = d2_new
!        psi_old  = psi_new
!
!     end if

  end do

  energy = energy/normalization
  accep  = dble(n_accep)/dble(nmax)

end subroutine pdmc



program qmc

  implicit none

  real(8)   , parameter :: a     = 1.2d0
  real(8)   , parameter :: dt    = 0.1d0
  real(8)   , parameter :: E_ref = -0.5d0
  real(8)   , parameter :: tau   = 100.d0
  integer(8), parameter :: nmax  = 100000
  integer   , parameter :: nruns = 30

  integer :: irun
  real(8) :: E(nruns), Acc(nruns)
  real(8) :: ave, err

  do irun = 1,nruns
     call pdmc(a, dt, nmax, E(irun), Acc(irun), tau, E_ref)
  enddo

  call ave_error(E,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err

  call ave_error(Acc,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err

end program qmc
