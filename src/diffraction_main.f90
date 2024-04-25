  module mod_Main   !tuy35 whole file
    use Size_FFT_Para
    implicit none

    integer trans

    real*8,parameter :: l0 = 1.d-9             !length unit
    real*8 os0   !tuy40
    character*8 :: passfilename

    !Order parameter and displacement
    real*8,allocatable,dimension(:,:,:,:),target :: oStruc   !tuy40 px,py,pz
    real*8,allocatable,dimension(:,:,:,:),target :: u
    real*8,target :: strainAvg(6)

    integer*8 kt0,ktMax                        !starting time step #; finishing time step #

    integer nPhase                             !total # of phases
    integer nStruc                             !total # of structural order parameters   !tuy40
    real*8,allocatable,dimension(:,:,:,:),target :: oPhase
!tuy39    real*8,allocatable,dimension(:,:,:) :: iPhase

    integer,allocatable,dimension(:) :: ixD,iyD,izD,outD,rankD
    real*8,allocatable,target :: IDiffr(:,:,:),DQ(:,:,:,:)
    real*8,target :: QCenter(3)
    ! Create a temp array DQ_add_Q save DQ + QCenter
    real*8,allocatable,target :: DQ_add_Q(:,:,:,:)
  end module



  program main

    ! use mod_interfaces
    ! use mod_fftw_mpi
    use mod_Main
    use mod_interface_diffraction
    use mod_mupro_base, only: type_mupro_SizeContext, mupro_size_setup
    use mod_mupro_fft, only: type_mupro_FFTContext, mupro_fft_setup
    use mod_mupro_io, only: mupro_output_4D_one_row, mupro_input_4D_one_row
    use mod_mupro_io, only: mupro_input_3D, mupro_output_3D

    implicit none
  
    integer i,j,k,l,m,n
    logical lexist
    type(type_mupro_SizeContext) sizeContext
    type(type_mupro_FFTContext) fftContext
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initiate mpi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call MPI_INIT(ierr)
    call MPI_comm_size(MPI_COMM_WORLD,process,ierr)
    call MPI_comm_rank(MPI_comm_world,rank,ierr)

4002  format("Hello, I am processor ",i3," of ",i2)
    do n = 0,process-1
      call MPI_Barrier(MPI_Comm_world,ierr)
      if(rank==n) then
!test        print 4002, rank,process-1
      endif
    enddo
    call sleep (1)
    call MPI_Barrier(MPI_Comm_world,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rank==0) then                                  
      open(unit = 1, file = "parameter.system.in")   !tuy38
      print *, "Input simulation system parameters"   !tuy40
      read(1,*)
      read(1,*),lx,ly,lz
      read(1,*),nx,ny,nz
      read(1,*),ns,nf
      read(1,*)
      read(1,*),nPhase   !tuy39
      read(1,*),nStruc   !tuy40
      read(1,*)
      read(1,*),strainAvg(1),strainAvg(2),strainAvg(3)
      read(1,*),strainAvg(4),strainAvg(5),strainAvg(6)
      close(1)
      print *   !tuy40
    endif

    call MPI_Barrier(MPI_Comm_world,ierr)

    call MPI_Bcast(lx,               1,   MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(ly,               1,   MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(lz,               1,   MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(nx,               1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(ny,               1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(nz,               1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(ns,               1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(nf,               1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(nPhase,           1,   MPI_integer, 0,MPI_Comm_World,ierr)   !tuy39
    call MPI_Bcast(nStruc,           1,   MPI_integer, 0,MPI_Comm_World,ierr)   !tuy40
    call MPI_Bcast(strainAvg,        6,   MPI_real8,   0,MPI_Comm_World,ierr)
 
    if (nPhase<1.or.nPhase>12) then   !tuy39b
      if (rank==0) print *, "This program allows 1~12 phases only."   !tuy40b
      if (rank==0) print *, "You cannot claim", nPhase, "phases."
      if (rank==0) print *, "Program cancelled"   !tuy40f
      goto 999
    endif   !tuy39f
    if (nStruc<0.or.nStruc>12) then   !tuy40b
      if (rank==0) print *, "This program allows 0~12 structural order parameters only."
      if (rank==0) print *, "You cannot claim", nStruc, "structural order parameters."
      if (rank==0) print *, "Program cancelled"
      goto 999
    endif   !tuy40f
    call MPI_Barrier(MPI_Comm_world,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! setup system size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dx=lx/real(nx)
    dy=ly/real(ny)
    dz=lz/real(nz)

    k1 = ns + 1
    k2 = ns + nf
    if (nf==0) then
      k1 = 1
      k2 = nz
    endif

    kt = kt0

	  ! call simSizeDeclare(nx,ny,nz,ns,nf,lx,ly,lz,kt_in=kt)
    sizeContext%nx = nx
    sizeContext%ny = ny
    sizeContext%nz = nz
    sizeContext%nf = nf
    sizecontext%ns = ns
    sizeContext%dx = dx
    sizeContext%dy = dy
    sizeContext%dz = dz
    call mupro_size_setup(sizeContext)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! call fourier_mpi_setup(R,C,HN,lstart,trans)
    call mupro_fft_setup(fftContext)

    Rn1 = fftContext%Rn1
    Rn2 = fftContext%Rn2
    Rn3 = fftContext%Rn3
    Cn1 = fftContext%Cn1
    Cn2 = fftContext%Cn2
    Cn3 = fftContext%Cn3
    Hn1 = fftContext%Hn1
    Hn2 = fftContext%Hn2
    lstartR = fftContext%lstart

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! allocate and initiate arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(oPhase(Rn3,Rn2,Rn1,nPhase));  oPhase=0.   !tuy39
!tuy39    allocate(iPhase(Rn3,Rn2,Rn1));         iPhase=0.

    allocate(u(Rn3,Rn2,Rn1,3));            u=0.

    allocate(oStruc(Rn3,Rn2,Rn1,nStruc));  oStruc=0.   !tuy40

!tuy40    allocate(px(Rn3,Rn2,Rn1));             px=0.
!    allocate(py(Rn3,Rn2,Rn1));             py=0.
!    allocate(pz(Rn3,Rn2,Rn1));             pz=0.

    allocate(IDiffr(Rn3,Rn2,Rn1));         IDiffr=0.
    allocate(DQ(Rn3,Rn2,Rn1,3));           DQ=0.
    allocate(DQ_add_Q(Rn3,Rn2,Rn1,3));     DQ_add_Q=0.

    call MPI_Barrier(Mpi_Comm_world,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input phase concentration field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lexist=.false.
    inquire(file='phaseFra.in',exist=lexist)

    if (lexist) then
      passfilename='phaseFra'
      call mupro_input_4D_one_row(passfilename, oPhase)
      ! if     (nPhase==1)  then; call InputxN_P(passfilename,oPhase(:,:,:,1))   !tuy39b
      ! elseif (nPhase==2)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2))
      ! elseif (nPhase==3)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3))
      ! elseif (nPhase==4)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4))
      ! elseif (nPhase==5)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5))
      ! elseif (nPhase==6)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6))
      ! elseif (nPhase==7)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6),oPhase(:,:,:,7))
      ! elseif (nPhase==8)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6),oPhase(:,:,:,7),oPhase(:,:,:,8))
      ! elseif (nPhase==9)  then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6),oPhase(:,:,:,7),oPhase(:,:,:,8),oPhase(:,:,:,9))
      ! elseif (nPhase==10) then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6),oPhase(:,:,:,7),oPhase(:,:,:,8),oPhase(:,:,:,9),oPhase(:,:,:,10))
      ! elseif (nPhase==11) then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6),oPhase(:,:,:,7),oPhase(:,:,:,8),oPhase(:,:,:,9),oPhase(:,:,:,10),oPhase(:,:,:,11))
      ! elseif (nPhase==12) then; call InputxN_P(passfilename,oPhase(:,:,:,1),oPhase(:,:,:,2),oPhase(:,:,:,3),oPhase(:,:,:,4),oPhase(:,:,:,5),oPhase(:,:,:,6),oPhase(:,:,:,7),oPhase(:,:,:,8),oPhase(:,:,:,9),oPhase(:,:,:,10),oPhase(:,:,:,11),oPhase(:,:,:,12))
      ! endif   !tuy39f

    else
      if(rank==0) print *, "File phaseFra.in not provided. Using a pure phase 1."   !tuy40
      if(rank==0) print *   !tuy40
      oPhase(:,:,:,1) = 1.d0   !tuy39

    endif

    oPhase(:k1-1,:,:,:)=0.
    oPhase(k2+1:,:,:,:)=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input structural order parameter field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    os0 = 1.d0   !tuy40
    lexist=.false.
    inquire(file='strucOrd.in',exist=lexist)

    if (lexist) then
      passfilename='strucOrd'
      if     (nStruc==0)  then; if(rank==0) print *, "File strucOrd.in skipped since nStruc = 0."; if(rank==0) print *   !tuy40b
      elseif  (nStruc/=0) then; call mupro_input_4D_one_row(passfilename, oStruc)
      ! elseif (nStruc==1)  then; call InputxN_P(passfilename,oStruc(:,:,:,1))
      ! elseif (nStruc==2)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2))
      ! elseif (nStruc==3)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3))
      ! elseif (nStruc==4)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4))
      ! elseif (nStruc==5)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5))
      ! elseif (nStruc==6)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6))
      ! elseif (nStruc==7)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6),oStruc(:,:,:,7))
      ! elseif (nStruc==8)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6),oStruc(:,:,:,7),oStruc(:,:,:,8))
      ! elseif (nStruc==9)  then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6),oStruc(:,:,:,7),oStruc(:,:,:,8),oStruc(:,:,:,9))
      ! elseif (nStruc==10) then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6),oStruc(:,:,:,7),oStruc(:,:,:,8),oStruc(:,:,:,9),oStruc(:,:,:,10))
      ! elseif (nStruc==11) then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6),oStruc(:,:,:,7),oStruc(:,:,:,8),oStruc(:,:,:,9),oStruc(:,:,:,10),oStruc(:,:,:,11))
      ! elseif (nStruc==12) then; call InputxN_P(passfilename,oStruc(:,:,:,1),oStruc(:,:,:,2),oStruc(:,:,:,3),oStruc(:,:,:,4),oStruc(:,:,:,5),oStruc(:,:,:,6),oStruc(:,:,:,7),oStruc(:,:,:,8),oStruc(:,:,:,9),oStruc(:,:,:,10),oStruc(:,:,:,11),oStruc(:,:,:,12))
      endif
!      call InputxN_P(passfilename,px,py,pz)
!      px = px /p0
!      py = py /p0
!      pz = pz /p0
      oStruc = oStruc /os0

    elseif (nStruc>0) then
      if (rank==0) print *, "Missing file strucOrd.in, which is required for nStruc > 0."
      if (rank==0) print *, "Program cancelled"
      goto 999
!      px = 0.d0
!      py = 0.d0
!tuy40f      pz = 0.d0

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input displacement field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lexist=.false.
    inquire(file='displace.in',exist=lexist)

    if (lexist) then
      passfilename='displace'
      call mupro_input_4D_one_row(passfilename, u)
      ! call InputxN_P(passfilename,u(:,:,:,1),u(:,:,:,2),u(:,:,:,3))

    else
      if(rank==0) print *, "File displace.in not provided. Using mechanical displacement = 0."   !tuy40
      if(rank==0) print *   !tuy40
      u = 0.d0

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! setup-routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call diffraction_setup(u,strainAvg,oPhase,oStruc,IDiffr,DQ,QCenter,os0,nPhase,nStruc,trans)   !tuy39   !tuy40
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! prepare for output !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!tuy37    if(rank==0) open(unit = 35, file = "qCenter.dat")
!    write(35,'("      kt    QCenter(1,2,3)(nm^-1)")')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kt = kt0
    call diffractionCalc

!tuy37    if(rank==0) write(35,1001) kt,QCenter(1),QCenter(2),QCenter(3)
1001  format(i10,30es16.7e3)

    passfilename = 'I'
    ! call outputxN_P(passfilename,IDiffr)
    call mupro_output_3D(passfilename, kt, IDiffr)

    passfilename = 'lg_{10}I'
    ! call outputxN_P(passfilename,dlog(IDiffr)/dlog(10.D0))
    call mupro_output_3D(passfilename, kt, dlog(IDiffr)/dlog(10.D0))

    passfilename = 'qVector'   !tuy37
    DQ_add_Q(:,:,:,1) = DQ(:,:,:,1) + QCenter(1)
    DQ_add_Q(:,:,:,2) = DQ(:,:,:,2) + QCenter(2)
    DQ_add_Q(:,:,:,3) = DQ(:,:,:,3) + QCenter(3)
    call mupro_output_4D_one_row(passfilename, kt, DQ_add_Q)
    ! call outputxN_P(passfilename,DQ(:,:,:,1)+QCenter(1),DQ(:,:,:,2)+QCenter(2),DQ(:,:,:,3)+QCenter(3))   !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_barrier(Mpi_comm_world,ierr)

!tuy37    if(rank==0) close(13)

	  print *, "Finished program on Rank ",rank

999 call MPI_finalize(ierr)

	end program
