subroutine variational_montecarlo(a,nmax,dt,energy,accep)

  implicit none

  real(8)   , intent(in)  :: a, dt
  integer(8), intent(in)  :: nmax 
  real(8)   , intent(out) :: energy, accep

  integer(8) :: istep
  integer(8) :: n_accep
  real(8)    :: sq_dt, chi(3), dq1_old, dq2_old, prod, u
  real(8)    :: psi_old, psi_new, dq1_new, dq2_new, argexpo, q
  real(8)    :: r1_old(3), r1_new(3), r2_old(3), r2_new(3)
  real(8)    :: d1_old(3), d1_new(3), d2_old(3), d2_new(3)

  double precision, external :: e_loc, psi

  sq_dt = dsqrt(dt)

  ! Initialization
  energy = 0.d0
  n_accep = 0_8

  call random_gauss(r1_old,3)
  call random_gauss(r2_old,3)

  call drift(a,r1_old,d1_old)
  call drift(a,r2_old,d2_old)
  dq1_old = d1_old(1)*d1_old(1) + &
            d1_old(2)*d1_old(2) + &
            d1_old(3)*d1_old(3)
  dq2_old = d2_old(1)*d2_old(1) + &
            d2_old(2)*d2_old(2) + &
            d2_old(3)*d2_old(3)

  psi_old = psi(a, r1_old, r2_old)

  do istep = 1,nmax

     energy = energy + e_loc(a,r1_old,r2_old)

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

     ! Metropolis
     prod = (d_new(1) + d_old(1))*(r1_new(1) - r1_old(1)) + &
            (d_new(2) + d_old(2))*(r1_new(2) - r1_old(2)) + &
            (d_new(3) + d_old(3))*(r1_new(3) - r1_old(3))

     argexpo = 0.5d0*(d2_new - d2_old)*dt + prod

     q = psi_new/psi_old
     q = dexp(-argexpo)*q*q

     call random_number(u)

     if (u.le.q) then

        n_accep = n_accep + 1_8

        r1_old(:) = r1_new(:)
        r2_old(:) = r2_new(:)
        d_old(:)  = d_new(:)
        d2_old    = d2_new
        psi_old   = psi_new

     end if

  end do

  energy = energy/dble(nmax)
  accep  = dble(n_accep)/dble(nmax)

end subroutine variational_montecarlo



program qmc

  implicit none

  real(8)   , parameter :: a     = 1.2d0
  real(8)   , parameter :: dt    = 1.d0
  integer(8), parameter :: nmax  = 100000
  integer   , parameter :: nruns = 30

  integer :: irun
  real(8) :: E(nruns), Acc(nruns)
  real(8) :: ave, err

  do irun = 1,nruns
    call variational_montecarlo(a,nmax,dt,E(irun),Acc(irun))
  enddo

  call ave_error(E,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err

  call ave_error(Acc,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err

end program qmc
