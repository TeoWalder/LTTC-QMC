!-----------------------------------------------!
!--------------- QMC Project -------------------!
!-----------------------------------------------!

program qmc

  implicit none
  ! Parameters
  real(8), parameter :: Ab    = 1.8897259886d0 ! Angstrom to Bohr conv
  ! Physical Variables
  integer              :: irun, ne, nn
  real(8)              :: ave, err
  real(8)              :: E_ref, a, Z = 1.d0
  real(8), allocatable :: Rn(:,:), E(:), Acc(:)
  real(8)              :: dt
  real(8)              :: tau
  ! Tecnical stuff
  integer(8)        :: nmax
  integer           :: nruns
  integer           :: i, ios
  character(len=3)  :: method, atom
  character(7)      :: inputfile = 'qmc.inp'
  character(len=10) :: filename

!----------------------------------------------------------!

!------ Initialization ------------------------------------!

  open(10, file=inputfile, status='old', iostat=ios)

    if (ios.ne.0) then
      write(*,*) '*** ERROR opening input file ***'
      stop
    else
      read(10,*) method                        ! var (vmc) or dif (pdmc)
      read(10,*) a                             ! Wavefunction Parameter
      read(10,*) ne                            ! n of Electrons
      read(10,*) filename                      ! Coordinate file (.xyz)
      read(10,*) dt                            ! Time Step
      read(10,*) nmax                          ! n of Steps
      read(10,*) nruns                         ! n of Walkers
      if (method.eq.'dif') then                ! only for pdmc:
        read(10,*) E_ref                       ! reference energy
        read(10,*) tau                         ! projection time 
      end if
    end if

  close(10)

  allocate(E(nruns), Acc(nruns))
  
  open(11, file=filename, status='old', iostat=ios)

    if (ios.ne.0) then
      write(*,*) '*** ERROR opening ', filename, '***'
      stop
    else
      read(11,*) nn
      read(11,*)

  allocate(Rn(nn,3))

      do i = 1,nn
        read(11,*) atom, Rn(i,:)
      end do
    end if

  close(11)

  Rn(:,:) = Ab*Rn(:,:)                   ! Conversion to a.u.

  if (atom.eq.'He') Z = 2.d0

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
    write(*,*) '*** ERROR, no such method ***'
    stop
  end if

!------ Output --------------------------------------------!

  call ave_error(E,nruns,ave,err)
  write(*,100) 'E = ', ave, ' +/- ', err

  call ave_error(Acc,nruns,ave,err)
  write(*,100) 'A = ', ave, ' +/- ', err

  open(12, file='e_loc.dat')
  do irun = 1,nruns
    write(12,*) irun, E(irun)
  end do
  close(12)

  100 format(a4,f15.9,a7,e15.6)

!----------------------------------------------------------!

  deallocate(Rn,E,Acc)
  stop

end program qmc
