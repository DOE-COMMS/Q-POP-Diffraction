  module mod_interface_diffraction   !tuy22 whole file

    implicit none
    #include"mpif.h"



    INTERFACE

      subroutine diffraction_setup(u,strainAvg,oPhase,oStruc,IDiffr,DQ,QCenter,os0_in,nPhase_in,nStruc_in,trans_in)   !tuy24   !tuy29   !tuy30   !tuy39   !tuy40
!tuy40        real*8,intent(IN),dimension(:,:,:),target :: px,py,pz
        real*8,intent(IN),target :: u(:,:,:,:),strainAvg(:)
        real*8,intent(IN),target :: oPhase(:,:,:,:)   !tuy29
        real*8,intent(IN),target :: oStruc(:,:,:,:)   !tuy40
        real*8,intent(IN),target :: IDiffr(:,:,:)   !tuy30
        real*8,intent(IN),target :: DQ(:,:,:,:),QCenter(:)   !tuy30
        real*8,intent(IN) :: os0_in   !tuy24   !tuy29
        integer,intent(IN) :: nPhase_in   !tuy39
        integer,intent(IN) :: nStruc_in   !tuy40
        integer,intent(IN) :: trans_in   !tuy24
      end subroutine

      subroutine ArrayFourierToRegular(ArrayA_in,ArrayB_in,Array_out,trans_in)
        real*8,intent(IN),dimension(:,:,:) :: ArrayA_in,ArrayB_in   !tuy23 ,allocatable
        real*8,intent(OUT),dimension(:,:,:) :: Array_out   !tuy30 ,allocatable
        integer,intent(IN) :: trans_in
      end subroutine

    END INTERFACE

  end module

  module Size_FFT_Para
    implicit none
    !system
    integer rank, process, ierr
    integer nx,ny,nz,nf,ns,k1,k2
    integer*8 kt
    real*8 lx,ly,lz,dx,dy,dz,dt0
    integer R(3), C(3), HN(2),lstart(3)
    integer Rn3,Rn2,Rn1
    integer Cn3,Cn2,Cn1
    integer Hn2,Hn1
    integer lstart3,lstart2,lstartR
    integer nn
    integer,allocatable,dimension(:) :: Rn1All,lstartRAll
  end module

  module mod_Diffraction   !tuy21   !tuy22moved

    implicit none

!tuy24    logical flagDiffraction                    !whether to output diffraction data
!tuy24    integer*8 ktOutDiffr                       !output step interval for diffraction data
!tuy33    integer hM1,hM2,hM3                        !reciprocal lattice point around which the diffraction pattern will be calculated (Miller index)   !tuy22
    real*8 q00(3)                                       !reciprocal lattice point around which the diffraction pattern will be calculated (rad.nm^-1)   !tuy33
    real*8 aC(3,3)                                      !lattice constants (nm)   !tuy32
    real*8,allocatable :: xAtom(:,:)                    !coordinates of all 5 atoms (unitless)   !tuy38b
    complex*16,allocatable :: fAtom(:,:)                !atomic form factors   !tuy29
    real*8,allocatable :: fAtomR(:,:),fAtomI(:,:)       !atomic form factors, real and imaginary parts   !tuy29
    real*8,allocatable,dimension(:,:,:) :: su1,su2,su3  !atomic-displacement-over-order-parameter coefficient (m^3/C)   !tuy39   !tuy40
!tuy38f    real*8 kScale                              !spatial scaling parameter (unitless)   !tuy24
  end module



  module diffraction

    implicit none

    real*8,parameter :: pi = dacos(-1.d0)
    real*8,parameter :: l0 = 1.d-9             !length unit   !tuy24
    integer nAtom                              !# of atoms in a unit cell
    integer nPhase                             !total # of phases   !tuy39
    integer nStruc                             !total # of structural order parameters   !tuy40
    integer trans
    real*8 bC(3,3)                             !reciprocal lattice constants (nm^-1)   !tuy32b
    real*8 VCell                               !unit cell volume
    real*8 DRAtom(5,3)                         !coordinates of all 5 atoms (nm)
    real*8 G10,G20,G30   !tuy32f
    integer hM1,hM2,hM3                        !"closet" Bragg peak to the calculated region (Miller index)   !tuy33
    real*8 G1,G2,G3   !tuy30   !tuy33
    real*8 qC10,qC20,qC30   !tuy33
    real*8,pointer :: qC1,qC2,qC3   !tuy33
    real*8,pointer,dimension(:,:,:) :: u1,u2,u3   !tuy40 p1,p2,p3,
    real*8,pointer,dimension(:) :: eAvg
    real*8,pointer,dimension(:,:,:,:) :: op   !tuy39  oPhase1,oPhase2   !tuy29
    real*8,pointer,dimension(:,:,:,:) :: os   !tuy40
    real*8,pointer,dimension(:,:,:) :: I_out   !tuy30
    real*8,pointer,dimension(:,:,:) :: dmqk1,dmqk2,dmqk3   !tuy30   !tuy33 All mqkG->mqk
    real*8,allocatable,dimension(:,:,:) :: mqkA1,mqkA2,mqkA3
    real*8,allocatable,dimension(:,:,:) :: mqkB1,mqkB2,mqkB3
    real*8,allocatable,dimension(:,:,:) :: mqk1_out,mqk2_out,mqk3_out   !tuy23 moved
    real*8,allocatable,dimension(:,:,:) :: region   !tuy33test
    complex*16,allocatable,dimension(:,:,:,:) :: fA,fB   !tuy23
    complex*16,allocatable,dimension(:,:,:) :: kfExpA,kfExpB,kfExp1,kfExp2,kfExp3,kTemp   !tuy41
    integer Cn1Sum
    integer,allocatable,dimension(:) :: Cn1All,lstartCAll

  end module



  subroutine diffraction_setup(u,strainAvg,oPhase,oStruc,IDiffr,DQ,QCenter,os0_in,nPhase_in,nStruc_in,trans_in)   !tuy24   !tuy29   !tuy30   !tuy39   !tuy40

    ! use mod_interfaces
    use Size_FFT_Para
    use mod_interface_diffraction
    ! use simSize
    ! use mod_fftw_mpi
    use mod_Diffraction
    use diffraction
    use mod_mupro_fft
    use mod_mupro_io, only: mupro_input_3D, mupro_output_3D

    implicit none

!tuy40    real*8,intent(IN),dimension(:,:,:),target :: px,py,pz
    real*8,intent(IN),target :: u(:,:,:,:),strainAvg(:)
    real*8,intent(IN),target :: oPhase(:,:,:,:)   !tuy29
    real*8,intent(IN),target :: oStruc(:,:,:,:)   !tuy40
    real*8,intent(IN),target :: IDiffr(:,:,:)   !tuy30
    real*8,intent(IN),target :: DQ(:,:,:,:),QCenter(:)   !tuy30
    real*8,intent(IN) :: os0_in   !tuy24   !tuy29
    integer,intent(IN) :: nPhase_in   !tuy39
    integer,intent(IN) :: nStruc_in   !tuy40
    integer,intent(IN) :: trans_in   !tuy24

    character*8 :: passfilename
    integer i,j,k,m,n   !tuy39
    logical lexist   !tuy35
    real*8 ix1,ix2,iy1,iy2,iz1,iz2   !tuy36
    real*8 tx,ty,tz   !tuy36
    real*8,allocatable,dimension(:) :: maxRegion   !tuy40

!tuy38    nAtom = 5   !tuy24
    nPhase = nPhase_in   !tuy39
    nStruc = nStruc_in   !tuy40
    trans = trans_in
!tuy40    p1 => px;           p2 => py;           p3 => pz
    u1 => u(:,:,:,1);   u2 => u(:,:,:,2);   u3 => u(:,:,:,3)
    eAvg => strainAvg
    op => oPhase   !tuy29   !tuy39
    os => oStruc   !tuy40
    I_out => IDiffr   !tuy30b
    dmqk1 => DQ(:,:,:,1);   dmqk2 => DQ(:,:,:,2);   dmqk3 => DQ(:,:,:,3)
    qC1 => QCenter(1);   qC2 => QCenter(2);   qC3 => QCenter(3)   !tuy33   !tuy30f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   !tuy24b partly moved        
    if (rank==0) then                 
      open(unit = 1, file = "parameter.atom.in")
      print *, "Input crystal information"   !tuy40
      read(1,*)
      read(1,*),q00(1),q00(2),q00(3)   !tuy33
      read(1,*),aC(1,1),aC(1,2),aC(1,3)   !tuy32
      read(1,*),aC(2,1),aC(2,2),aC(2,3)   !
      read(1,*),aC(3,1),aC(3,2),aC(3,3)   !
      read(1,*),nAtom   !tuy38b
    endif

    call MPI_Barrier(MPI_Comm_world,ierr)

    call MPI_Bcast(q00,     3,   MPI_real8,   0,MPI_Comm_World,ierr)   !tuy33
    call MPI_Bcast(aC,      9,   MPI_real8,   0,MPI_Comm_World,ierr)   !tuy32
    call MPI_Bcast(nAtom,   1,   MPI_integer, 0,MPI_Comm_World,ierr)

    allocate(xAtom(nAtom,3))
    allocate(fAtom(nPhase,nAtom))   !tuy39b
    allocate(fAtomR(nPhase,nAtom),fAtomI(nPhase,nAtom))
    allocate(su1(nStruc,nPhase,nAtom),su2(nStruc,nPhase,nAtom),su3(nStruc,nPhase,nAtom))   !tuy40
!tuy40    allocate(su11(nPhase,nAtom),su12(nPhase,nAtom),su13(nPhase,nAtom))   !tuy40
!    allocate(su21(nPhase,nAtom),su22(nPhase,nAtom),su23(nPhase,nAtom))
!    allocate(su31(nPhase,nAtom),su32(nPhase,nAtom),su33(nPhase,nAtom))   !tuy39f

    call MPI_Barrier(MPI_Comm_world,ierr)

    if (rank==0) then                 
      print *, "Input atom coordinates"   !tuy40
      read(1,*)
      do n = 1, nAtom   !tuy39 All i->n for Atom ID
        read(1,*),xAtom(n,1),xAtom(n,2),xAtom(n,3)
      enddo

      do i = 1,nPhase   !tuy39b
        read(1,*)
        read(1,*),m
        print *, "Input parameters for phase #", m   !tuy40
        do n = 1, nAtom
          read(1,*),fAtomR(m,n),fAtomI(m,n)
        enddo
        do n = 1, nAtom
          do k = 1, nStruc   !tuy40b
            read(1,*),su1(k,m,n),su2(k,m,n),su3(k,m,n)
!          read(1,*),su11(m,n),su12(m,n),su13(m,n)
!          read(1,*),su21(m,n),su22(m,n),su23(m,n)
!          read(1,*),su31(m,n),su32(m,n),su33(m,n)
          enddo   !tuy40f
        enddo
      enddo   !tuy39f

!      read(1,*),kScale
      close(1)
      print *   !tuy40
    endif

    call MPI_Barrier(MPI_Comm_world,ierr)

    call MPI_Bcast(xAtom,  nAtom*3,             MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(fAtomR, nAtom*nPhase,        MPI_real8,   0,MPI_Comm_World,ierr)   !tuy29   !tuy39b
    call MPI_Bcast(fAtomI, nAtom*nPhase,        MPI_real8,   0,MPI_Comm_World,ierr)   !
    call MPI_Bcast(su1,    nAtom*nPhase*nStruc, MPI_real8,   0,MPI_Comm_World,ierr)   !tuy40b
    call MPI_Bcast(su2,    nAtom*nPhase*nStruc, MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(su3,    nAtom*nPhase*nStruc, MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su11,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su12,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su13,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su21,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su22,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su23,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su31,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!    call MPI_Bcast(su32,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)
!tuy40f    call MPI_Bcast(su33,    nAtom*nPhase,MPI_real8,   0,MPI_Comm_World,ierr)   !tuy39f
!tuy38f    call MPI_Bcast(kScale,  1,   MPI_real8,   0,MPI_Comm_World,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! unitless constants and variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!tuy38    fAtomR(:,4) = fAtomR(:,3)   !tuy29
!    fAtomI(:,4) = fAtomI(:,3)   !
!    fAtomR(:,5) = fAtomR(:,3)   !
!    fAtomI(:,5) = fAtomI(:,3)   !
    fAtom = fAtomR + (0,1)*fAtomI
!    if (rank==0) then   !test
!      write(*,*), "Atomic form factors"
!      write(*,*), fAtom(1,1)
!      write(*,*), fAtom(1,2)
!      write(*,*), fAtom(1,3)
!      write(*,*), fAtom(1,4)
!      write(*,*), fAtom(1,5)
!      write(*,*), fAtom(2,1)
!      write(*,*), fAtom(2,2)
!      write(*,*), fAtom(2,3)
!      write(*,*), fAtom(2,4)
!      write(*,*), fAtom(2,5)
!    endif

!tuy38b    su22(1) = su11(1);  su33(1) = su11(1)
!    su22(2) = su11(2);  su33(2) = su11(2)
!    su22(3) = su11(4);  su33(3) = su11(4)
!    su22(4) = su11(3);  su33(4) = su11(4)
!    su22(5) = su11(4);  su33(5) = su11(3)

!    su11(5) = su11(4)

    su1 = su1 / (l0/os0_in)   !tuy40b
    su2 = su2 / (l0/os0_in)
    su3 = su3 / (l0/os0_in)
!    su11 = su11 / (l0/p0_in)
!    su12 = su12 / (l0/p0_in)
!    su13 = su13 / (l0/p0_in)
!    su21 = su21 / (l0/p0_in)
!    su22 = su22 / (l0/p0_in)
!    su23 = su23 / (l0/p0_in)
!    su31 = su31 / (l0/p0_in)
!    su32 = su32 / (l0/p0_in)
!tuy40b    su33 = su33 / (l0/p0_in)   !tuy38f

    call MPI_Barrier(MPI_Comm_world,ierr)

    if(.not.allocated(region)) allocate(region(Rn3,Rn2,Rn1))   !tuy33testb

    lexist=.false.   !tuy35b
    inquire(file='region.in',exist=lexist)

    if (lexist) then
      passfilename='region'
      call mupro_input_3D(passfilename, region)
      ! call InputxN_P(passfilename,region)   !tuy33testf


    else
      if(rank==0) print *, "File region.in not provided. Using default region."   !tuy40
      if(rank==0) print *   !tuy40b
      ix1 = (nx+1)/2.d0 - nx/pi - lstartR   !tuy36b
      ix2 = (nx+1)/2.d0 + nx/pi - lstartR
      iy1 = (ny+1)/2.d0 - ny/pi
      iy2 = (ny+1)/2.d0 + ny/pi
      iz1 = (nz+1)/2.d0 - nz/pi
      iz2 = (nz+1)/2.d0 + nz/pi   !tuy40f
      tx = nx/8.d0
      ty = ny/8.d0
      tz = nz/8.d0
      
      region = 1.d0
      if(nf==0) then; do i = 1,Rn3;  region(i,:,:) =                 ( tanh((i-iz1)/tz) - tanh((i-iz2)/tz) );  enddo;   endif
                      do i = 1,Rn2;  region(:,i,:) = region(:,i,:) * ( tanh((i-iy1)/ty) - tanh((i-iy2)/ty) );  enddo
                      do i = 1,Rn1;  region(:,:,i) = region(:,:,i) * ( tanh((i-ix1)/tx) - tanh((i-ix2)/tx) );  enddo

      if(.not.allocated(maxRegion)) allocate(maxRegion(0:process-1))   !tuy40b
      maxRegion(rank) = maxval(region)
      call MPI_Barrier(MPI_Comm_world,ierr)

      do n = 0,process-1
        call MPI_Bcast(maxRegion(n),     1,MPI_real8  ,n,MPI_Comm_World,ierr)
      enddo
      call MPI_Barrier(MPI_Comm_world,ierr)

      region = region/maxval(maxRegion)   !tuy40f

      passfilename='region'
      ! call outputxN_P(passfilename,region)   !tuy36f
      call mupro_output_3D(passfilename, kt, region)

    endif   !tuy35f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   !tuy24f
    if(.not.allocated(mqkA1)) allocate(mqkA1(Cn3,Cn2,Cn1));         mqkA1=0.
    if(.not.allocated(mqkA2)) allocate(mqkA2(Cn3,Cn2,Cn1));         mqkA2=0.
    if(.not.allocated(mqkA3)) allocate(mqkA3(Cn3,Cn2,Cn1));         mqkA3=0.
    if(.not.allocated(mqkB1)) allocate(mqkB1(Cn3,Cn2,Cn1));         mqkB1=0.
    if(.not.allocated(mqkB2)) allocate(mqkB2(Cn3,Cn2,Cn1));         mqkB2=0.
    if(.not.allocated(mqkB3)) allocate(mqkB3(Cn3,Cn2,Cn1));         mqkB3=0.
    if(.not.allocated(fA))    allocate(fA(Cn3,Cn2,Cn1,nAtom));      fA=(0.,0.)   !tuy23
    if(.not.allocated(fB))    allocate(fB(Cn3,Cn2,Cn1,nAtom));      fB=(0.,0.)   !
    if(.not.allocated(kfExpA))allocate(kfExpA(Cn3,Cn2,Cn1));        kfExpA=(0.,0.)   !tuy41b
    if(.not.allocated(kfExpB))allocate(kfExpB(Cn3,Cn2,Cn1));        kfExpB=(0.,0.)
    if(.not.allocated(kfExp1))allocate(kfExp1(Cn3,Cn2,Cn1));        kfExp1=(0.,0.)
    if(.not.allocated(kfExp2))allocate(kfExp2(Cn3,Cn2,Cn1));        kfExp2=(0.,0.)
    if(.not.allocated(kfExp3))allocate(kfExp3(Cn3,Cn2,Cn1));        kfExp3=(0.,0.)
    if(.not.allocated(kTemp)) allocate(kTemp(Cn3,Cn2,Cn1));         kTemp=(0.,0.)   !tuy41f

    VCell = aC(1,1)*aC(2,2)*aC(3,3) + aC(1,2)*aC(2,3)*aC(3,1) + aC(1,3)*aC(2,1)*aC(3,2)&   !tuy32b
          - aC(1,1)*aC(2,3)*aC(3,2) - aC(1,2)*aC(2,1)*aC(3,3) - aC(1,3)*aC(2,2)*aC(3,1)

    bc(1,1) = 2*pi/VCell * (aC(2,2)*aC(3,3)-aC(2,3)*aC(3,2))
    bc(1,2) = 2*pi/VCell * (aC(2,3)*aC(3,1)-aC(2,1)*aC(3,3))
    bc(1,3) = 2*pi/VCell * (aC(2,1)*aC(3,2)-aC(2,2)*aC(3,1))

    bc(2,1) = 2*pi/VCell * (aC(1,3)*aC(3,2)-aC(1,2)*aC(3,3))
    bc(2,2) = 2*pi/VCell * (aC(1,1)*aC(3,3)-aC(1,3)*aC(3,1))
    bc(2,3) = 2*pi/VCell * (aC(1,2)*aC(3,1)-aC(1,1)*aC(3,2))

    bc(3,1) = 2*pi/VCell * (aC(1,2)*aC(2,3)-aC(1,3)*aC(2,2))
    bc(3,2) = 2*pi/VCell * (aC(1,3)*aC(2,1)-aC(1,1)*aC(2,3))
    bc(3,3) = 2*pi/VCell * (aC(1,1)*aC(2,2)-aC(1,2)*aC(2,1))

    hM1 = idnint((q00(1)*aC(1,1) + q00(2)*aC(1,2) + q00(3)*aC(1,3))/(2*pi))   !tuy33
    hM2 = idnint((q00(1)*aC(2,1) + q00(2)*aC(2,2) + q00(3)*aC(2,3))/(2*pi))   !
    hM3 = idnint((q00(1)*aC(3,1) + q00(2)*aC(3,2) + q00(3)*aC(3,3))/(2*pi))   !

    G10 = hM1*bc(1,1) + hM2*bc(2,1) + hM3*bc(3,1)
    G20 = hM1*bc(1,2) + hM2*bc(2,2) + hM3*bc(3,2)
    G30 = hM1*bc(1,3) + hM2*bc(2,3) + hM3*bc(3,3)

    qC10 = G10!tuy42 + idnint((q00(1)-G10)/(2*pi/lx)) * 2*pi/lx   !tuy33b
    qC20 = G20! + idnint((q00(2)-G20)/(2*pi/ly)) * 2*pi/ly
    qC30 = G30! + idnint((q00(3)-G30)/(2*pi/lz)) * 2*pi/lz

!    G1 = G10
!    G2 = G20
!    G3 = G30   !tuy32f
!    if (rank==0) print *, "G1 =", G10  !test   !tuy40b
!    if (rank==0) print *, "G2 =", G20
!    if (rank==0) print *, "G3 =", G30
!    if (rank==0) print *   !tuy40f

    mqkA1 =  mqk1_3 + qC10   !Here, mqkA & mqkB is q
    mqkA2 =  mqk2_3 + qC20
    mqkA3 =  mqk3_3 + qC30

    mqkB1 = -mqk1_3 + qC10
    mqkB2 = -mqk2_3 + qC20
    mqkB3 = -mqk3_3 + qC30

    if(mod(nx,2)==0) then
      if(trans==1) then;  i=nx/2+1;                               mqkB1(:,i,:) = mqk1_3(:,i,:) + qC10
      else;               i=nx/2+1-lstart3;  if (i>=0.and.i<=Cn1) mqkB1(:,:,i) = mqk1_3(:,:,i) + qC10;  endif
    endif

    if(mod(ny,2)==0) then
      if(trans==1) then;  j=ny/2+1-lstart3;  if (j>=0.and.j<=Cn1) mqkB2(:,:,j) = mqk2_3(:,:,j) + qC20
      else;               j=ny/2+1;                               mqkB2(:,j,:) = mqk2_3(:,j,:) + qC20;  endif
    endif

    if(mod(nz,2)==0) then
      k=nz/2+1;   mqkB3(k,:,:) = mqk3_3(k,:,:) + qC30   !tuy33f
    endif

    if(.not.allocated(mqk1_out)) allocate(mqk1_out(Rn3,Rn2,Rn1));  mqk1_out=0.
    if(.not.allocated(mqk2_out)) allocate(mqk2_out(Rn3,Rn2,Rn1));  mqk2_out=0.
    if(.not.allocated(mqk3_out)) allocate(mqk3_out(Rn3,Rn2,Rn1));  mqk3_out=0.

!test    write(*,*) "trans =",trans

    call ArrayFourierToRegular(mqkA1,mqkB1,mqk1_out,trans)
    call ArrayFourierToRegular(mqkA2,mqkB2,mqk2_out,trans)
    call ArrayFourierToRegular(mqkA3,mqkB3,mqk3_out,trans)

!tuy30    passfilename = 'mqk0'
!    call outputxN_P(passfilename,mqk1_out,mqk2_out,mqk3_out)
!tuy23    deallocate(mqk1_out, mqk2_out, mqk3_out)

    kTemp=(0,1)*(mqkA1-G10)*dx;  kfExp1=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp1=1.d0   !tuy41b
    kTemp=(0,1)*(mqkA2-G20)*dy;  kfExp2=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp2=1.d0
    kTemp=(0,1)*(mqkA3-G30)*dz;  kfExp3=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp3=1.d0
    kfExpA = kfExp1*kfExp2*kfExp3   !tuy41f

    kTemp=(0,1)*(mqkB1-G10)*dx;  kfExp1=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp1=1.d0
    kTemp=(0,1)*(mqkB2-G20)*dy;  kfExp2=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp2=1.d0
    kTemp=(0,1)*(mqkB3-G30)*dz;  kfExp3=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp3=1.d0
    kfExpB = kfExp1*kfExp2*kfExp3   !tuy41f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate structural factor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   !tuy24 moved
    do n=1,nAtom   !tuy23b
      DRAtom(n,1) = xAtom(n,1)*ac(1,1) + xAtom(n,2)*ac(2,1) + xAtom(n,3)*ac(3,1)   !tuy32b
      DRAtom(n,2) = xAtom(n,1)*ac(1,2) + xAtom(n,2)*ac(2,2) + xAtom(n,3)*ac(3,2)
      DRAtom(n,3) = xAtom(n,1)*ac(1,3) + xAtom(n,2)*ac(2,3) + xAtom(n,3)*ac(3,3)
      fA(:,:,:,n) = 1./VCell * cdexp( -(0,1) * (mqkA1*DRAtom(n,1) + mqkA2*DRAtom(n,2) + mqkA3*DRAtom(n,3)) )   !tuy29
      fB(:,:,:,n) = 1./VCell * cdexp( -(0,1) * (mqkB1*DRAtom(n,1) + mqkB2*DRAtom(n,2) + mqkB3*DRAtom(n,3)) )   !   !tuy32f
    enddo   !tuy23f

    mqkA1 = mqkA1 - qC10   !tuy33b   !tuy27b   !From this moment on, mqkA & mqkB becomes /Delta q
    mqkA2 = mqkA2 - qC20
    mqkA3 = mqkA3 - qC30

    mqkB1 = mqkB1 - qC10
    mqkB2 = mqkB2 - qC20
    mqkB3 = mqkB3 - qC30

    dmqk1 = mqk1_out - qC10   !tuy30b
    dmqk2 = mqk2_out - qC20
    dmqk3 = mqk3_out - qC30   !tuy33f   !dmqk is also /Delta q

!tuy37    passfilename = 'mqkD'
!    call outputxN_P(passfilename,dmqk1,dmqk2,dmqk3)   !tuy27f   !tuy30f

  end subroutine



  subroutine diffractionCalc

    ! use mod_interfaces
    use mod_interface_diffraction
    ! use simSize
    ! use mod_fftw_mpi
    use mod_Diffraction
    use diffraction
    use Size_FFT_Para
    use mod_mupro_fft, only: mupro_fft_forward

    implicit none

    real*8,allocatable,dimension(:,:,:) :: tempR   !tuy27b
    complex*16,allocatable,dimension(:,:,:) :: tempC
!    complex*16,allocatable,dimension(:,:,:) :: pk1,pk2,pk3   !tuy23b partly moved
!    complex*16,allocatable,dimension(:,:,:) :: uk1,uk2,uk3
    real*8,allocatable,dimension(:,:,:,:) :: uNM1,uNM2,uNM3   !tuy39 All uN->uNM
    real*8,allocatable,dimension(:,:,:) :: q0Gr   !tuy33
    complex*16,allocatable,dimension(:,:,:,:) :: fExpNM  !tuy39
    complex*16,allocatable,dimension(:,:,:,:) :: fExpN   !tuy29b
    complex*16,allocatable,dimension(:,:,:,:) :: fU1ExpN,fU2ExpN,fU3ExpN
    complex*16,allocatable,dimension(:,:,:,:) :: fExpNRk,fExpNIk
    complex*16,allocatable,dimension(:,:,:,:) :: fU1ExpNRk,fU2ExpNRk,fU3ExpNRk
    complex*16,allocatable,dimension(:,:,:,:) :: fU1ExpNIk,fU2ExpNIk,fU3ExpNIk   !tuy29f
    complex*16,allocatable,dimension(:,:,:) :: AmpA,AmpB
    real*8,allocatable,dimension(:,:,:) :: IA,IB   !tuy30,I_out
    character*8 :: passfilename   !tuy23f
    integer i,j,k,l,m,n

    if(.not.allocated(tempR))    allocate(tempR(Rn3,Rn2,Rn1));          tempR=0.   !tuy27b
    if(.not.allocated(tempC))    allocate(tempC(Cn3,Cn2,Cn1));          tempC=0.
    if(.not.allocated(uNM1))     allocate(uNM1(Rn3,Rn2,Rn1,nAtom));     uNM1=(0.,0.)
    if(.not.allocated(uNM2))     allocate(uNM2(Rn3,Rn2,Rn1,nAtom));     uNM2=(0.,0.)
    if(.not.allocated(uNM3))     allocate(uNM3(Rn3,Rn2,Rn1,nAtom));     uNM3=(0.,0.)
    if(.not.allocated(q0Gr))     allocate(q0Gr(Rn3,Rn2,Rn1));           q0Gr=0.   !tuy33
    if(.not.allocated(fExpNM))   allocate(fExpNM(Rn3,Rn2,Rn1,nAtom));   fExpNM=(0.,0.)   !tuy39
    if(.not.allocated(fExpN))    allocate(fExpN(Rn3,Rn2,Rn1,nAtom));    fExpN=(0.,0.)   !tuy29b
    if(.not.allocated(fU1ExpN))  allocate(fU1ExpN(Rn3,Rn2,Rn1,nAtom));  fU1ExpN=(0.,0.)
    if(.not.allocated(fU2ExpN))  allocate(fU2ExpN(Rn3,Rn2,Rn1,nAtom));  fU2ExpN=(0.,0.)
    if(.not.allocated(fU3ExpN))  allocate(fU3ExpN(Rn3,Rn2,Rn1,nAtom));  fU3ExpN=(0.,0.)   !tuy29f
!    if(.not.allocated(pk1))      allocate(pk1(Cn3,Cn2,Cn1));            pk1=(0.,0.)   !tuy23b partly moved
!    if(.not.allocated(pk2))      allocate(pk2(Cn3,Cn2,Cn1));            pk2=(0.,0.)
!    if(.not.allocated(pk3))      allocate(pk3(Cn3,Cn2,Cn1));            pk3=(0.,0.)
!    if(.not.allocated(uk1))      allocate(uk1(Cn3,Cn2,Cn1));            uk1=(0.,0.)
!    if(.not.allocated(uk2))      allocate(uk2(Cn3,Cn2,Cn1));            uk2=(0.,0.)
!    if(.not.allocated(uk3))      allocate(uk3(Cn3,Cn2,Cn1));            uk3=(0.,0.)
    if(.not.allocated(fExpNRk))  allocate(fExpNRk(Cn3,Cn2,Cn1,nAtom));  fExpNRk=(0.,0.)   !tuy29b
    if(.not.allocated(fU1ExpNRk))allocate(fU1ExpNRk(Cn3,Cn2,Cn1,nAtom));fU1ExpNRk=(0.,0.)
    if(.not.allocated(fU2ExpNRk))allocate(fU2ExpNRk(Cn3,Cn2,Cn1,nAtom));fU2ExpNRk=(0.,0.)
    if(.not.allocated(fU3ExpNRk))allocate(fU3ExpNRk(Cn3,Cn2,Cn1,nAtom));fU3ExpNRk=(0.,0.)
    if(.not.allocated(fExpNIk))  allocate(fExpNIk(Cn3,Cn2,Cn1,nAtom));  fExpNIk=(0.,0.)
    if(.not.allocated(fU1ExpNIk))allocate(fU1ExpNIk(Cn3,Cn2,Cn1,nAtom));fU1ExpNIk=(0.,0.)
    if(.not.allocated(fU2ExpNIk))allocate(fU2ExpNIk(Cn3,Cn2,Cn1,nAtom));fU2ExpNIk=(0.,0.)
    if(.not.allocated(fU3ExpNIk))allocate(fU3ExpNIk(Cn3,Cn2,Cn1,nAtom));fU3ExpNIk=(0.,0.)   !tuy29f
!tuy27f
    if(.not.allocated(AmpA))     allocate(AmpA(Cn3,Cn2,Cn1));           AmpA=(0.,0.)
    if(.not.allocated(AmpB))     allocate(AmpB(Cn3,Cn2,Cn1));           AmpB=(0.,0.)
    if(.not.allocated(IA))       allocate(IA(Cn3,Cn2,Cn1));             IA=0.
    if(.not.allocated(IB))       allocate(IB(Cn3,Cn2,Cn1));             IB=0.
!tuy30    if(.not.allocated(I_out))    allocate(I_out(Rn3,Rn2,Rn1));          I_out=0.   !tuy23f

    call MPI_Bcast(eAvg,3,MPI_real8,0,MPI_Comm_World,ierr)   !tuy25


!tuy32b    G1 = hM1 * 2*pi / (aC(1)*(1+eAvg(1)))   !tuy23b !For nonzero average shear strains, need to do a lattice vector
!    G2 = hM2 * 2*pi / (aC(2)*(1+eAvg(2)))
!    G3 = hM3 * 2*pi / (aC(3)*(1+eAvg(3)))
    G1 = G10 - qC10*eAvg(1) - qC20*eAvg(6) - qC30*eAvg(5)   !tuy33b
    G2 = G20 - qC10*eAvg(6) - qC20*eAvg(2) - qC30*eAvg(4)
    G3 = G30 - qC10*eAvg(5) - qC20*eAvg(4) - qC30*eAvg(3)   !tuy32f
!test    if(rank==0) write(*,*) "eAvg =",eAvg(1),eAvg(2),eAvg(3)

    qC1 = qC10 + G1 - G10
    qC2 = qC20 + G2 - G20
    qC3 = qC30 + G3 - G30

    mqk1_out = dmqk1 + qC1   !tuy30b
    mqk2_out = dmqk2 + qC2
    mqk3_out = dmqk3 + qC3

!    passfilename = 'mqk'
!    call outputxN_P(passfilename,mqk1_out,mqk2_out,mqk3_out)   !tuy23f   !tuy30f

    q0Gr = 0.
    do i=1,Rn1;  q0Gr(:,:,i)=q0Gr(:,:,i)+(qC1-G1)*(i+lstartR)*dx;  enddo   !tuy40
    do i=1,Rn2;  q0Gr(:,i,:)=q0Gr(:,i,:)+(qC2-G2)*i*dy;  enddo
    do i=1,Rn3;  q0Gr(i,:,:)=q0Gr(i,:,:)+(qC3-G3)*i*dz;  enddo   !tuy33f

!tuy27b
!    call forward_mpi(p1,pk1)
!    call forward_mpi(p2,pk2)
!    call forward_mpi(p3,pk3)
!    call forward_mpi(u1,uk1)
!    call forward_mpi(u2,uk2)
!    call forward_mpi(u3,uk3)

    fExpN = 0.d0   !tuy39b
    fU1ExpN = 0.d0
    fU2ExpN = 0.d0
    fU3ExpN = 0.d0

    do n = 1,nAtom
      do m = 1,nPhase
        uNM1(:,:,:,n) = u1(:,:,:)   !tuy40b
        uNM2(:,:,:,n) = u2(:,:,:)
        uNM3(:,:,:,n) = u3(:,:,:)
        do k = 1,nStruc
          uNM1(:,:,:,n) = uNM1(:,:,:,n) + su1(k,m,n)*os(:,:,:,k)
          uNM2(:,:,:,n) = uNM2(:,:,:,n) + su2(k,m,n)*os(:,:,:,k)
          uNM3(:,:,:,n) = uNM3(:,:,:,n) + su3(k,m,n)*os(:,:,:,k)
!          uNM1(:,:,:,n) = uNM1(:,:,:,n) + su11(m,n)*p1 + su12(m,n)*p2 + su13(m,n)*p3   !tuy38
!          uNM2(:,:,:,n) = uNM2(:,:,:,n) + su21(m,n)*p1 + su22(m,n)*p2 + su23(m,n)*p3   !
!          uNM3(:,:,:,n) = uNM3(:,:,:,n) + su31(m,n)*p1 + su32(m,n)*p2 + su33(m,n)*p3   !
        enddo   !tuy40f
        fExpNM(:,:,:,n) = cdexp( -(0,1)* (qC1*uNM1(:,:,:,n) + qC2*uNM2(:,:,:,n) + qC3*uNM3(:,:,:,n) + q0Gr) ) *fAtom(m,n)*op(:,:,:,m)*region   !tuy41
          fExpN(:,:,:,n) =   fExpN(:,:,:,n) + fExpNM(:,:,:,n)
        fU1ExpN(:,:,:,n) = fU1ExpN(:,:,:,n) + uNM1(:,:,:,n) * fExpNM(:,:,:,n)
        fU2ExpN(:,:,:,n) = fU2ExpN(:,:,:,n) + uNM2(:,:,:,n) * fExpNM(:,:,:,n)
        fU3ExpN(:,:,:,n) = fU3ExpN(:,:,:,n) + uNM3(:,:,:,n) * fExpNM(:,:,:,n)
      enddo   !tuy39f
    enddo

    do n = 1,nAtom
      tempR = real(  fExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);    fExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(  fExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);    fExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = real(fU1ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU1ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(fU1ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU1ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = real(fU2ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU2ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(fU2ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU2ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = real(fU3ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU3ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(fU3ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU3ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)   !tuy29f
    enddo

!tuy27f

!test    n=1
!    call ArrayFourierToRegular(real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)), real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quReal1'
!    call outputxN_P(passfilename,I_out)
!    call ArrayFourierToRegular(imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),-imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quImag1'
!    call outputxN_P(passfilename,I_out)
!
!    n=2
!    call ArrayFourierToRegular(real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)), real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quReal2'
!    call outputxN_P(passfilename,I_out)
!    call ArrayFourierToRegular(imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),-imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quImag2'
!    call outputxN_P(passfilename,I_out)
!
!    n=3
!    call ArrayFourierToRegular(real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)), real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quReal3'
!    call outputxN_P(passfilename,I_out)
!    call ArrayFourierToRegular(imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),-imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quImag3'
!    call outputxN_P(passfilename,I_out)
!
!    n=4
!    call ArrayFourierToRegular(real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)), real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quReal4'
!    call outputxN_P(passfilename,I_out)
!    call ArrayFourierToRegular(imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),-imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quImag4'
!    call outputxN_P(passfilename,I_out)
!
!    n=5
!    call ArrayFourierToRegular(real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)), real(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quReal5'
!    call outputxN_P(passfilename,I_out)
!    call ArrayFourierToRegular(imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),-imag(mqkA1*uNk1(:,:,:,n) + mqkA2*uNk2(:,:,:,n) + mqkA3*uNk3(:,:,:,n)),I_out,trans)
!    passfilename = 'quImag5'
!    call outputxN_P(passfilename,I_out)

    AmpA = (0.,0.)   !tuy23b
    AmpB = (0.,0.)
    do n=1,nAtom
      AmpA = AmpA + fA(:,:,:,n) * (                          fExpNRk(:,:,:,n)  + (0,1)*          fExpNIk(:,:,:,n)   &   !tuy27b   !tuy29b
                                    -(0,1)* (mqkA1*       fU1ExpNRk(:,:,:,n)  + mqkA2*       fU2ExpNRk(:,:,:,n)  + mqkA3*       fU3ExpNRk(:,:,:,n)) &
                                    +       (mqkA1*       fU1ExpNIk(:,:,:,n)  + mqkA2*       fU2ExpNIk(:,:,:,n)  + mqkA3*       fU3ExpNIk(:,:,:,n))  )
      AmpB = AmpB + fB(:,:,:,n) * (                 dconjg(  fExpNRk(:,:,:,n)) + (0,1)* dconjg(  fExpNIk(:,:,:,n))  &
                                    -(0,1)* (mqkB1*dconjg(fU1ExpNRk(:,:,:,n)) + mqkB2*dconjg(fU2ExpNRk(:,:,:,n)) + mqkB3*dconjg(fU3ExpNRk(:,:,:,n)))&
                                    +       (mqkB1*dconjg(fU1ExpNIk(:,:,:,n)) + mqkB2*dconjg(fU2ExpNIk(:,:,:,n)) + mqkB3*dconjg(fU3ExpNIk(:,:,:,n))) )   !tuy29f
    enddo

!    if(rank==0) then   !tuyissue singular point
!      do n=1,nAtom
!        AmpA(1,1,1) = AmpA(1,1,1) + fA(1,1,1,n)   !tuy24  !tuy26
!        AmpB(1,1,1) = AmpB(1,1,1) + fB(1,1,1,n)   !       !
!      enddo
!tuy27f    endif

    AmpA = AmpA * kfExpA   !tuy41
    AmpB = AmpB * kfExpB   !

    IA = max(cdabs(AmpA)**2,1.d-50)
    IB = max(cdabs(AmpB)**2,1.d-50)
!tuy38    if(rank==0) then   !tuy27 tuyissue singular point
!      IA(1,1,1) = IA(1,1,1) * kScale
!      IB(1,1,1) = IB(1,1,1) * kScale
!    endif   !tuy27f

    call ArrayFourierToRegular(IA,IB,I_out,trans)
!tuy30    passfilename = 'diffrInt'
!    call outputxN_P(passfilename,I_out)

!tuy27b    deallocate(pk1,pk2,pk3,uk1,uk2,uk3,uNk1,uNk2,uNk3)
    deallocate(tempR,tempC)
    deallocate(q0Gr)   !tuy33
    deallocate(uNM1,uNM2,uNM3,fExpNM)   !tuy29b   !tuy39
    deallocate(fExpN,fU1ExpN,fU2ExpN,fU3ExpN)   !tuy39
    deallocate(fExpNRk,fU1ExpNRk,fU2ExpNRk,fU3ExpNRk)
    deallocate(fExpNIk,fU1ExpNIk,fU2ExpNIk,fU3ExpNIk)   !tuy27f   !tuy29f
    deallocate(AmpA,AmpB,IA,IB)   !tuy23f   !tuy30 ,I_out

  end subroutine



  subroutine ArrayFourierToRegular(ArrayA_in,ArrayB_in,Array_out,trans_in)

    ! use simSize
    ! use mod_fftw_mpi
    use mod_interface_diffraction
    use Size_FFT_Para

    implicit none

    real*8,intent(IN),dimension(:,:,:) :: ArrayA_in,ArrayB_in   !tuy23 ,allocatable
    real*8,intent(OUT),dimension(:,:,:) :: Array_out   !tuy30 ,allocatable
    integer,intent(IN) :: trans_in

    integer Cn1Sum
    integer,allocatable,dimension(:) :: Cn1All,lstartCAll
    real*8,allocatable,dimension(:,:,:) :: Array_All,ArrayA_All,ArrayB_All,Array_out_All
    integer i,j,k,l,m,n,ii,jj,kk

    if(.not.allocated(Cn1All))     allocate(Cn1All(0:process-1));      Cn1All=0.
    if(.not.allocated(lstartCAll)) allocate(lstartCAll(0:process-1));  lstartCAll=0.
    Cn1All(rank) = Cn1
    lstartCAll(rank) = lstart3
    call MPI_Barrier(MPI_Comm_world,ierr)

    do n = 0,process-1
      call MPI_Bcast(Cn1All(n),     1,MPI_integer,n,MPI_Comm_World,ierr)
      call MPI_Bcast(lstartCAll(n), 1,MPI_integer,n,MPI_Comm_World,ierr)
    enddo
    call MPI_Barrier(MPI_Comm_world,ierr)

    Cn1Sum = sum(Cn1All)
    call MPI_Barrier(MPI_Comm_world,ierr)

567  format ("rank =",i3,",   Cn1Sum =",i3,",   Cn1 =",i3,",   Cn2 =",i3,",   Cn3 =",i3,",   lstartC =",i3,",   lstart3 =",i3)
!    write(*,567),rank,Cn1Sum,Cn1,Cn2,Cn3,lstartCAll(rank),lstart3

    if(.not.allocated(ArrayA_All))    allocate(ArrayA_All(Cn3,Cn2,Cn1Sum));  ArrayA_All=0.
    if(.not.allocated(ArrayB_All))    allocate(ArrayB_All(Cn3,Cn2,Cn1Sum));  ArrayB_All=0.
    if(.not.allocated(Array_All))     allocate(Array_All(nz,ny,nx));         Array_All=0.
    if(.not.allocated(Array_out_All)) allocate(Array_out_All(nz,ny,nx));     Array_out_All=0.
!tuy30    if(.not.allocated(Array_out))     allocate(Array_out(Rn3,Rn2,Rn1));      Array_out=0.

    call MPI_Allgatherv(ArrayA_in(:,:,:),Cn1*Cn2*Cn3,MPI_real8,ArrayA_All(:,:,:),Cn1All(:)*Cn2*Cn3,lstartCAll(:)*Cn2*Cn3,MPI_real8,MPI_Comm_World,ierr)
    call MPI_Allgatherv(ArrayB_in(:,:,:),Cn1*Cn2*Cn3,MPI_real8,ArrayB_All(:,:,:),Cn1All(:)*Cn2*Cn3,lstartCAll(:)*Cn2*Cn3,MPI_real8,MPI_Comm_World,ierr)

    do i = 1,nx
    do j = 1,ny
    do k = 1,nz/2+1
      kk = mod(nz+2-k-1+nz,nz)+1
      jj = mod(ny+2-j-1+ny,ny)+1
      ii = mod(nx+2-i-1+nx,nx)+1
      if (trans_in==1) then
        Array_All(k,j,i) = ArrayA_All(k,i,j)
        if(kk>nz/2+1 .and. kk<=nz .and. jj<=ny .and. ii<=nx) Array_All(kk,jj,ii) = ArrayB_All(k,i,j)
      else
        Array_All(k,j,i) = ArrayA_All(k,j,i)
        if(kk>nz/2+1 .and. kk<=nz .and. jj<=ny .and. ii<=nx) Array_All(kk,jj,ii) = ArrayB_All(k,j,i)
      endif
    enddo
    enddo
    enddo

!    write(*,*) "Array_All size =", size(Array_All,1), size(Array_All,2), size(Array_All,3)
    do i = 1,nx
    do j = 1,ny
    do k = 1,nz
      ii = mod(i+nx/2-1,nx)+1
      jj = mod(j+ny/2-1,ny)+1
      kk = mod(k+nz/2-1,nz)+1
      Array_out_All(kk,jj,ii) = Array_All(k,j,i)
    enddo
    enddo
    enddo

    do i=1,Rn1
      Array_out(:,:,i) = Array_out_All(:,:,i+lstartR)
    enddo

    deallocate(ArrayA_All, ArrayB_All, Array_All, Array_out_All)

  end subroutine
