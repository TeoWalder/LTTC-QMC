program qmc

  implicit none

  real(8)     , parameter :: dt      = 0.1d0
  real(8)     , parameter :: tau     = 100.d0
  integer(8)  , parameter :: nmax    = 100000
  integer     , parameter :: nruns   = 30
  character(5), parameter :: inpfile = "input"

  integer              :: irun, ne, nn
  real(8)              :: E(nruns), Acc(nruns)
  real(8)              :: ave, err
  real(8)              :: E_ref, a, Z = 1.d0
  real(8), allocatable :: Rn(:,:)

  integer           :: i, ios
  character(len=3)  :: method, atom
  character(len=10) :: filename

  open(10, file=inpfile, status="old", iostat=ios)

    if (ios.ne.0) then
      write(*,*) "*** ERROR opening input file ***"
      stop
    else
      read(10,*) method
      read(10,*) a
      read(10,*) ne
      read(10,*) filename
      read(10,*) E_ref
    end if

  close(10)
  
  open(11, file=filename, status="old", iostat=ios)

    if (ios.ne.0) then
      write(*,*) "*** ERROR opening xyz file ***"
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

  if (atom.eq."He") Z = 2.d0

  if (method.eq."var") then
    do irun = 1,nruns
      call VMC(a, nmax, dt, E(irun), Acc(irun), Rn, ne, nn, Z)
    enddo
  else if (method.eq."dif") then
    do irun = 1,nruns
       call PDMC(a, dt, nmax, E(irun), Acc(irun), tau, E_ref, Rn, ne, nn, Z)
    enddo
  else
    write(*,*) "*** ERROR, no such method ***"
    stop
  end if

  call ave_error(E,nruns,ave,err)
  write(*,100) "E = ", ave, " +/- ", err

  call ave_error(Acc,nruns,ave,err)
  write(*,100) "A = ", ave, " +/- ", err

  100 format(a4,f15.9,a5,e15.6)

end program qmc
