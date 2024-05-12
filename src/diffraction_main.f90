  module mod_Main
    use Size_FFT_Para
    implicit none

    integer trans

    real*8,parameter :: l0 = 1.d-9             !length unit
    real*8 os0

    character*8 :: passfilename

    !Order parameter and displacement
    real*8,allocatable,dimension(:,:,:,:),target :: oStruc
    real*8,allocatable,dimension(:,:,:,:),target :: u
    real*8,target :: strainAvg(6)

    integer*8 kt0,ktMax                        !starting time step #; finishing time step #

    integer nPhase                             !total # of phases
    integer nStruc                             !total # of structural order parameters
    real*8,allocatable,dimension(:,:,:,:),target :: oPhase

    integer,allocatable,dimension(:) :: ixD,iyD,izD,outD,rankD
    real*8,allocatable,target :: IDiffr(:,:,:),DQ(:,:,:,:)
    real*8,target :: QCenter(3)
    ! Create a temp array DQ_add_Q save DQ + QCenter
    real*8,allocatable,target :: DQ_add_Q(:,:,:,:)
  end module



  program main

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
      open(unit = 1, file = "parameter.system.in")
      print *, "Input simulation system parameters"
      read(1,*)
      read(1,*),lx,ly,lz
      read(1,*),nx,ny,nz
      read(1,*),ns,nf
      read(1,*)
      read(1,*),nPhase
      read(1,*),nStruc
      read(1,*)
      read(1,*),strainAvg(1),strainAvg(2),strainAvg(3)
      read(1,*),strainAvg(4),strainAvg(5),strainAvg(6)
      close(1)
      print *
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
    call MPI_Bcast(nPhase,           1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(nStruc,           1,   MPI_integer, 0,MPI_Comm_World,ierr)
    call MPI_Bcast(strainAvg,        6,   MPI_real8,   0,MPI_Comm_World,ierr)
 
    if (nPhase<1.or.nPhase>12) then
      if (rank==0) print *, "This program allows 1~12 phases only."
      if (rank==0) print *, "You cannot claim", nPhase, "phases."
      if (rank==0) print *, "Program cancelled"
      goto 999
    endif
    if (nStruc<0.or.nStruc>12) then
      if (rank==0) print *, "This program allows 0~12 structural order parameters only."
      if (rank==0) print *, "You cannot claim", nStruc, "structural order parameters."
      if (rank==0) print *, "Program cancelled"
      goto 999
    endif
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
    allocate(oPhase(nPhase,Rn3,Rn2,Rn1));  oPhase=0.

    allocate(u(3,Rn3,Rn2,Rn1));            u=0.

    allocate(oStruc(nStruc,Rn3,Rn2,Rn1));  oStruc=0.

    allocate(IDiffr(Rn3,Rn2,Rn1));         IDiffr=0.
    allocate(DQ(3,Rn3,Rn2,Rn1));           DQ=0.
    allocate(DQ_add_Q(3,Rn3,Rn2,Rn1));     DQ_add_Q=0.

    call MPI_Barrier(Mpi_Comm_world,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input phase concentration field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lexist=.false.
    inquire(file='phaseFra.in',exist=lexist)

    if (lexist) then
      passfilename='phaseFra'
      call mupro_input_4D_one_row(passfilename, oPhase)

    else
      if(rank==0) print *, "File phaseFra.in not provided. Using a pure phase 1."
      if(rank==0) print *
      oPhase(1,:,:,:) = 1.d0

    endif

    oPhase(:,:k1-1,:,:)=0.
    oPhase(:,k2+1:,:,:)=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input structural order parameter field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    os0 = 1.d0
    lexist=.false.
    inquire(file='strucOrd.in',exist=lexist)

    if (lexist) then
      passfilename='strucOrd'
      if     (nStruc==0)  then; if(rank==0) print *, "File strucOrd.in skipped since nStruc = 0."; if(rank==0) print *
      elseif  (nStruc/=0) then; call mupro_input_4D_one_row(passfilename, oStruc)
      endif
      oStruc = oStruc /os0

    elseif (nStruc>0) then
      if (rank==0) print *, "Missing file strucOrd.in, which is required for nStruc > 0."
      if (rank==0) print *, "Program cancelled"
      goto 999

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input displacement field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lexist=.false.
    inquire(file='displace.in',exist=lexist)

    if (lexist) then
      passfilename='displace'
      call mupro_input_4D_one_row(passfilename, u)

    else
      if(rank==0) print *, "File displace.in not provided. Using mechanical displacement = 0."
      if(rank==0) print *
      u = 0.d0

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! setup-routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call diffraction_setup(u,strainAvg,oPhase,oStruc,IDiffr,DQ,QCenter,os0,nPhase,nStruc,trans)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kt = kt0
    call diffractionCalc

    passfilename = 'I'
    ! call outputxN_P(passfilename,IDiffr)
    call mupro_output_3D(passfilename, kt, IDiffr)

    passfilename = 'lg_{10}I'
    ! call outputxN_P(passfilename,dlog(IDiffr)/dlog(10.D0))
    call mupro_output_3D(passfilename, kt, dlog(IDiffr)/dlog(10.D0))

    passfilename = 'qVector'
    DQ_add_Q(1,:,:,:) = DQ(1,:,:,:) + QCenter(1)
    DQ_add_Q(2,:,:,:) = DQ(2,:,:,:) + QCenter(2)
    DQ_add_Q(3,:,:,:) = DQ(3,:,:,:) + QCenter(3)
    call mupro_output_4D_one_row(passfilename, kt, DQ_add_Q)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_barrier(Mpi_comm_world,ierr)

	  print *, "Finished program on Rank ",rank

999 call MPI_finalize(ierr)

	end program
