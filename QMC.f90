!-----------------------------------------------!
!--------------- QMC Project -------------------!
!-----------------------------------------------!

program qmc

  implicit none
  ! Parameters
  real(8), parameter :: Ab = 1.8897259886d0  ! Angstrom to Bohr conversion
  ! Physical variables
  integer              :: irun, ne, nn
  real(8)              :: ave, err
  real(8)              :: E_ref, a
  real(8), allocatable :: Rn(:,:), E(:), Acc(:), Z(:)
  real(8)              :: dt
  real(8)              :: tau
  real(8)              :: CC(3)
  ! Tecnical variables
  integer(8)        :: nmax
  integer           :: nruns
  integer           :: i, ios
  character(len=3)  :: method, atom
  character(7)      :: inputfile = 'qmc.inp'
  character(len=10) :: geom

!----------------------------------------------------------!

!------ Initialization ------------------------------------!

  open(10, file=inputfile, status='old', iostat=ios)

    if (ios.ne.0) stop '*** ERROR opening input file'

    read(10,*) method                    ! var (vmc) or dif (pdmc)
    read(10,*) a                         ! Wavefunction Parameter
    read(10,*) ne                        ! n of Electrons
    read(10,*) geom                      ! Coordinate file (.xyz)
    read(10,*) dt                        ! Time Step
    read(10,*) nmax                      ! n of Steps
    read(10,*) nruns                     ! n of Walkers
    if (method.eq.'dif') then
      read(10,*)
      read(10,*)
      read(10,*) E_ref                   ! reference energy
      read(10,*) tau                     ! projection time 
    end if

  close(10)

  allocate(E(nruns), Acc(nruns))

  ! read geometry

  open(11, file=geom, status='old', iostat=ios)

    if (ios.ne.0) stop '*** ERROR opening xyz file'

    read(11,*) nn
    read(11,*)

  allocate(Rn(nn,3),Z(nn))

    do i = 1,nn
      read(11,*) atom, Rn(i,:)

      if (atom.eq.'H') then
        Z(i) = 1.d0
      else if (atom.eq.'He') then 
        Z(i) = 2.d0
      else 
        stop 'Invalid atom'
      end if

    end do

  close(11)

  ! a.u. conversion
  Rn(:,:) = Ab*Rn(:,:)

  ! compute Charge Center
  CC(:) = 0.d0
  do i = 1,nn
    CC(:) = CC(:) + Z(i)*Rn(i,:)
  end do
  CC(:) = CC(:)/sum(Z(:))

  ! translate the geometry by the mass center
  Rn(:,1) = Rn(:,1) - CC(1)
  Rn(:,2) = Rn(:,2) - CC(2)
  Rn(:,3) = Rn(:,3) - CC(3)


!------ Monte Carlo Calculation ---------------------------!

  if (method.eq.'var') then
    do irun = 1,nruns
      call VMC(a, nmax, dt, E(irun), Acc(irun), Rn, ne, nn, Z)
    enddo
  else if (method.eq.'dif') then
    do irun = 1,nruns
      call PDMC(a, dt, nmax, E(irun), Acc(irun), tau, E_ref, Rn, ne, nn, Z)
    enddo
  else
    stop '*** ERROR, no such method'
  end if


!------ Output --------------------------------------------!

  call ave_error(E,nruns,ave,err)
  write(*,100) 'E = ', ave, ' +/- ', err

  call ave_error(Acc,nruns,ave,err)
  write(*,100) 'A = ', ave, ' +/- ', err

  ! Plot local energy
  open(12, file='e_loc.dat')
    do irun = 1,nruns
      write(12,*) irun, E(irun)
    end do
  close(12)

  100 format(a4,f15.9,a7,e15.6)

!----------------------------------------------------------!

  deallocate(Rn,E,Acc,Z)

end program qmc
