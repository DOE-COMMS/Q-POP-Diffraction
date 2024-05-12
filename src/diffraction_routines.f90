  module mod_interface_diffraction

    implicit none
    !include"mpif.h"



    INTERFACE

      subroutine diffraction_setup(u,strainAvg,oPhase,oStruc,IDiffr,DQ,QCenter,os0_in,nPhase_in,nStruc_in,trans_in)
        real*8,intent(IN),target :: u(:,:,:,:),strainAvg(:)
        real*8,intent(IN),target :: oPhase(:,:,:,:)
        real*8,intent(IN),target :: oStruc(:,:,:,:)
        real*8,intent(IN),target :: IDiffr(:,:,:)
        real*8,intent(IN),target :: DQ(:,:,:,:),QCenter(:)
        real*8,intent(IN) :: os0_in
        integer,intent(IN) :: nPhase_in
        integer,intent(IN) :: nStruc_in
        integer,intent(IN) :: trans_in
      end subroutine

      subroutine ArrayFourierToRegular(ArrayA_in,ArrayB_in,Array_out,trans_in)
        real*8,intent(IN),dimension(:,:,:) :: ArrayA_in,ArrayB_in
        real*8,intent(OUT),dimension(:,:,:) :: Array_out
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

  module mod_Diffraction

    implicit none

    real*8 q00(3)                                       !reciprocal lattice point around which the diffraction pattern will be calculated (rad.nm^-1)
    real*8 aC(3,3)                                      !lattice constants (nm)
    real*8,allocatable :: xAtom(:,:)                    !coordinates of all n atoms (unitless)
    complex*16,allocatable :: fAtom(:,:)                !atomic form factors
    real*8,allocatable :: fAtomR(:,:),fAtomI(:,:)       !atomic form factors, real and imaginary parts
    real*8,allocatable,dimension(:,:,:) :: su1,su2,su3  !atomic-displacement-over-order-parameter coefficient (m^3/C)
  end module



  module diffraction

    implicit none

    real*8,parameter :: pi = dacos(-1.d0)
    real*8,parameter :: l0 = 1.d-9             !length unit
    integer nAtom                              !# of atoms in a unit cell
    integer nPhase                             !total # of phases
    integer nStruc                             !total # of structural order parameters
    integer trans
    real*8 bC(3,3)                             !reciprocal lattice constants (nm^-1)
    real*8 VCell                               !unit cell volume
    real*8,allocatable :: DRAtom(:,:)          !coordinates of all n atoms (nm)
    real*8 G10,G20,G30
    integer hM1,hM2,hM3                        !"closet" Bragg peak to the calculated region (Miller index)
    real*8 G1,G2,G3
    real*8 qC10,qC20,qC30
    real*8,pointer :: qC1,qC2,qC3
    real*8,pointer,dimension(:,:,:) :: u1,u2,u3
    real*8,pointer,dimension(:) :: eAvg
    real*8,pointer,dimension(:,:,:,:) :: op
    real*8,pointer,dimension(:,:,:,:) :: os
    real*8,pointer,dimension(:,:,:) :: I_out
    real*8,pointer,dimension(:,:,:) :: dmqk1,dmqk2,dmqk3
    real*8,allocatable,dimension(:,:,:) :: mqkA1,mqkA2,mqkA3
    real*8,allocatable,dimension(:,:,:) :: mqkB1,mqkB2,mqkB3
    real*8,allocatable,dimension(:,:,:) :: mqk1_out,mqk2_out,mqk3_out
    real*8,allocatable,dimension(:,:,:) :: region
    complex*16,allocatable,dimension(:,:,:,:) :: fA,fB
    complex*16,allocatable,dimension(:,:,:) :: kfExpA,kfExpB,kfExp1,kfExp2,kfExp3,kTemp
    integer Cn1Sum
    integer,allocatable,dimension(:) :: Cn1All,lstartCAll

  end module



  subroutine diffraction_setup(u,strainAvg,oPhase,oStruc,IDiffr,DQ,QCenter,os0_in,nPhase_in,nStruc_in,trans_in)

    use Size_FFT_Para
    use mod_interface_diffraction
    use mod_Diffraction
    use diffraction
    use mod_mupro_fft
    use mod_mupro_io, only: mupro_input_3D, mupro_output_3D

    implicit none

    real*8,intent(IN),target :: u(:,:,:,:),strainAvg(:)
    real*8,intent(IN),target :: oPhase(:,:,:,:)
    real*8,intent(IN),target :: oStruc(:,:,:,:)
    real*8,intent(IN),target :: IDiffr(:,:,:)
    real*8,intent(IN),target :: DQ(:,:,:,:),QCenter(:)
    real*8,intent(IN) :: os0_in
    integer,intent(IN) :: nPhase_in
    integer,intent(IN) :: nStruc_in
    integer,intent(IN) :: trans_in

    character*8 :: passfilename
    integer i,j,k,m,n
    logical lexist
    real*8 ix1,ix2,iy1,iy2,iz1,iz2
    real*8 tx,ty,tz
    real*8,allocatable,dimension(:) :: maxRegion

    nPhase = nPhase_in
    nStruc = nStruc_in
    trans = trans_in
    u1 => u(1,:,:,:);   u2 => u(2,:,:,:);   u3 => u(3,:,:,:)
    eAvg => strainAvg
    op => oPhase
    os => oStruc
    I_out => IDiffr
    dmqk1 => DQ(1,:,:,:);   dmqk2 => DQ(2,:,:,:);   dmqk3 => DQ(3,:,:,:)
    qC1 => QCenter(1);   qC2 => QCenter(2);   qC3 => QCenter(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! input parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rank==0) then
      open(unit = 1, file = "parameter.atom.in")
      print *, "Input crystal information"
      read(1,*)
      read(1,*),q00(1),q00(2),q00(3)
      read(1,*),aC(1,1),aC(1,2),aC(1,3)
      read(1,*),aC(2,1),aC(2,2),aC(2,3)
      read(1,*),aC(3,1),aC(3,2),aC(3,3)
      read(1,*),nAtom
    endif

    call MPI_Barrier(MPI_Comm_world,ierr)

    call MPI_Bcast(q00,     3,   MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(aC,      9,   MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(nAtom,   1,   MPI_integer, 0,MPI_Comm_World,ierr)

    allocate(xAtom(nAtom,3))
    allocate(fAtom(nPhase,nAtom))
    allocate(fAtomR(nPhase,nAtom),fAtomI(nPhase,nAtom))
    allocate(su1(nStruc,nPhase,nAtom),su2(nStruc,nPhase,nAtom),su3(nStruc,nPhase,nAtom))

    call MPI_Barrier(MPI_Comm_world,ierr)

    if (rank==0) then                 
      print *, "Input atom coordinates"
      read(1,*)
      do n = 1, nAtom
        read(1,*),xAtom(n,1),xAtom(n,2),xAtom(n,3)
      enddo

      do i = 1,nPhase
        read(1,*)
        read(1,*),m
        print *, "Input parameters for phase #", m
        do n = 1, nAtom
          read(1,*),fAtomR(m,n),fAtomI(m,n)
        enddo
        do n = 1, nAtom
          do k = 1, nStruc
            read(1,*),su1(k,m,n),su2(k,m,n),su3(k,m,n)
          enddo
        enddo
      enddo

      close(1)
      print *
    endif

    call MPI_Barrier(MPI_Comm_world,ierr)

    call MPI_Bcast(xAtom,  nAtom*3,             MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(fAtomR, nAtom*nPhase,        MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(fAtomI, nAtom*nPhase,        MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(su1,    nAtom*nPhase*nStruc, MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(su2,    nAtom*nPhase*nStruc, MPI_real8,   0,MPI_Comm_World,ierr)
    call MPI_Bcast(su3,    nAtom*nPhase*nStruc, MPI_real8,   0,MPI_Comm_World,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! unitless constants and variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fAtom = fAtomR + (0,1)*fAtomI
    su1 = su1 / (l0/os0_in)
    su2 = su2 / (l0/os0_in)
    su3 = su3 / (l0/os0_in)

    call MPI_Barrier(MPI_Comm_world,ierr)

    if(.not.allocated(region)) allocate(region(Rn3,Rn2,Rn1))

    lexist=.false.
    inquire(file='region.in',exist=lexist)

    if (lexist) then
      passfilename='region'
      call mupro_input_3D(passfilename, region)

    else
      if(rank==0) print *, "File region.in not provided. Using default region."
      if(rank==0) print *
      ix1 = (nx+1)/2.d0 - nx/pi - lstartR
      ix2 = (nx+1)/2.d0 + nx/pi - lstartR
      iy1 = (ny+1)/2.d0 - ny/pi
      iy2 = (ny+1)/2.d0 + ny/pi
      iz1 = (nz+1)/2.d0 - nz/pi
      iz2 = (nz+1)/2.d0 + nz/pi
      tx = nx/8.d0
      ty = ny/8.d0
      tz = nz/8.d0
      
      region = 1.d0
      if(nf==0) then; do i = 1,Rn3;  region(i,:,:) =                 ( tanh((i-iz1)/tz) - tanh((i-iz2)/tz) );  enddo;   endif
                      do i = 1,Rn2;  region(:,i,:) = region(:,i,:) * ( tanh((i-iy1)/ty) - tanh((i-iy2)/ty) );  enddo
                      do i = 1,Rn1;  region(:,:,i) = region(:,:,i) * ( tanh((i-ix1)/tx) - tanh((i-ix2)/tx) );  enddo

      if(.not.allocated(maxRegion)) allocate(maxRegion(0:process-1))
      maxRegion(rank) = maxval(region)
      call MPI_Barrier(MPI_Comm_world,ierr)

      do n = 0,process-1
        call MPI_Bcast(maxRegion(n),     1,MPI_real8  ,n,MPI_Comm_World,ierr)
      enddo
      call MPI_Barrier(MPI_Comm_world,ierr)

      region = region/maxval(maxRegion)

      passfilename='region'
      call mupro_output_3D(passfilename, kt, region)

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(.not.allocated(mqkA1)) allocate(mqkA1(Cn3,Cn2,Cn1));         mqkA1=0.
    if(.not.allocated(mqkA2)) allocate(mqkA2(Cn3,Cn2,Cn1));         mqkA2=0.
    if(.not.allocated(mqkA3)) allocate(mqkA3(Cn3,Cn2,Cn1));         mqkA3=0.
    if(.not.allocated(mqkB1)) allocate(mqkB1(Cn3,Cn2,Cn1));         mqkB1=0.
    if(.not.allocated(mqkB2)) allocate(mqkB2(Cn3,Cn2,Cn1));         mqkB2=0.
    if(.not.allocated(mqkB3)) allocate(mqkB3(Cn3,Cn2,Cn1));         mqkB3=0.
    if(.not.allocated(fA))    allocate(fA(Cn3,Cn2,Cn1,nAtom));      fA=(0.,0.)
    if(.not.allocated(fB))    allocate(fB(Cn3,Cn2,Cn1,nAtom));      fB=(0.,0.)
    if(.not.allocated(kfExpA))allocate(kfExpA(Cn3,Cn2,Cn1));        kfExpA=(0.,0.)
    if(.not.allocated(kfExpB))allocate(kfExpB(Cn3,Cn2,Cn1));        kfExpB=(0.,0.)
    if(.not.allocated(kfExp1))allocate(kfExp1(Cn3,Cn2,Cn1));        kfExp1=(0.,0.)
    if(.not.allocated(kfExp2))allocate(kfExp2(Cn3,Cn2,Cn1));        kfExp2=(0.,0.)
    if(.not.allocated(kfExp3))allocate(kfExp3(Cn3,Cn2,Cn1));        kfExp3=(0.,0.)
    if(.not.allocated(kTemp)) allocate(kTemp(Cn3,Cn2,Cn1));         kTemp=(0.,0.)

    VCell = aC(1,1)*aC(2,2)*aC(3,3) + aC(1,2)*aC(2,3)*aC(3,1) + aC(1,3)*aC(2,1)*aC(3,2)&
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

    hM1 = idnint((q00(1)*aC(1,1) + q00(2)*aC(1,2) + q00(3)*aC(1,3))/(2*pi))
    hM2 = idnint((q00(1)*aC(2,1) + q00(2)*aC(2,2) + q00(3)*aC(2,3))/(2*pi))
    hM3 = idnint((q00(1)*aC(3,1) + q00(2)*aC(3,2) + q00(3)*aC(3,3))/(2*pi))

    G10 = hM1*bc(1,1) + hM2*bc(2,1) + hM3*bc(3,1)
    G20 = hM1*bc(1,2) + hM2*bc(2,2) + hM3*bc(3,2)
    G30 = hM1*bc(1,3) + hM2*bc(2,3) + hM3*bc(3,3)

    qC10 = G10
    qC20 = G20
    qC30 = G30

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
      k=nz/2+1;   mqkB3(k,:,:) = mqk3_3(k,:,:) + qC30
    endif

    if(.not.allocated(mqk1_out)) allocate(mqk1_out(Rn3,Rn2,Rn1));  mqk1_out=0.
    if(.not.allocated(mqk2_out)) allocate(mqk2_out(Rn3,Rn2,Rn1));  mqk2_out=0.
    if(.not.allocated(mqk3_out)) allocate(mqk3_out(Rn3,Rn2,Rn1));  mqk3_out=0.

!test    write(*,*) "trans =",trans

    call ArrayFourierToRegular(mqkA1,mqkB1,mqk1_out,trans)
    call ArrayFourierToRegular(mqkA2,mqkB2,mqk2_out,trans)
    call ArrayFourierToRegular(mqkA3,mqkB3,mqk3_out,trans)

    kTemp=(0,1)*(mqkA1-G10)*dx;  kfExp1=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp1=1.d0
    kTemp=(0,1)*(mqkA2-G20)*dy;  kfExp2=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp2=1.d0
    kTemp=(0,1)*(mqkA3-G30)*dz;  kfExp3=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp3=1.d0
    kfExpA = kfExp1*kfExp2*kfExp3

    kTemp=(0,1)*(mqkB1-G10)*dx;  kfExp1=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp1=1.d0
    kTemp=(0,1)*(mqkB2-G20)*dy;  kfExp2=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp2=1.d0
    kTemp=(0,1)*(mqkB3-G30)*dz;  kfExp3=1.d0/kTemp*(cdexp(kTemp/2.d0)-cdexp(-kTemp/2.d0));  where(cdabs(kTemp)<=1.d-8) kfExp3=1.d0
    kfExpB = kfExp1*kfExp2*kfExp3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate structural factor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(.not.allocated(DRAtom)) allocate(DRAtom(nAtom,3));  DRAtom=0.

    do n=1,nAtom
      DRAtom(n,1) = xAtom(n,1)*ac(1,1) + xAtom(n,2)*ac(2,1) + xAtom(n,3)*ac(3,1)
      DRAtom(n,2) = xAtom(n,1)*ac(1,2) + xAtom(n,2)*ac(2,2) + xAtom(n,3)*ac(3,2)
      DRAtom(n,3) = xAtom(n,1)*ac(1,3) + xAtom(n,2)*ac(2,3) + xAtom(n,3)*ac(3,3)
      fA(:,:,:,n) = 1./VCell * cdexp( -(0,1) * (mqkA1*DRAtom(n,1) + mqkA2*DRAtom(n,2) + mqkA3*DRAtom(n,3)) )
      fB(:,:,:,n) = 1./VCell * cdexp( -(0,1) * (mqkB1*DRAtom(n,1) + mqkB2*DRAtom(n,2) + mqkB3*DRAtom(n,3)) )
    enddo

    mqkA1 = mqkA1 - qC10    !From this moment on, mqkA & mqkB becomes /Delta q
    mqkA2 = mqkA2 - qC20
    mqkA3 = mqkA3 - qC30

    mqkB1 = mqkB1 - qC10
    mqkB2 = mqkB2 - qC20
    mqkB3 = mqkB3 - qC30

    dmqk1 = mqk1_out - qC10
    dmqk2 = mqk2_out - qC20
    dmqk3 = mqk3_out - qC30   !dmqk is also /Delta q

  end subroutine



  subroutine diffractionCalc

    use mod_interface_diffraction
    use mod_Diffraction
    use diffraction
    use Size_FFT_Para
    use mod_mupro_fft, only: mupro_fft_forward

    implicit none

    real*8,allocatable,dimension(:,:,:) :: tempR
    complex*16,allocatable,dimension(:,:,:) :: tempC
    real*8,allocatable,dimension(:,:,:,:) :: uNM1,uNM2,uNM3
    real*8,allocatable,dimension(:,:,:) :: q0Gr
    complex*16,allocatable,dimension(:,:,:,:) :: fExpNM
    complex*16,allocatable,dimension(:,:,:,:) :: fExpN
    complex*16,allocatable,dimension(:,:,:,:) :: fU1ExpN,fU2ExpN,fU3ExpN
    complex*16,allocatable,dimension(:,:,:,:) :: fExpNRk,fExpNIk
    complex*16,allocatable,dimension(:,:,:,:) :: fU1ExpNRk,fU2ExpNRk,fU3ExpNRk
    complex*16,allocatable,dimension(:,:,:,:) :: fU1ExpNIk,fU2ExpNIk,fU3ExpNIk
    complex*16,allocatable,dimension(:,:,:) :: AmpA,AmpB
    real*8,allocatable,dimension(:,:,:) :: IA,IB
    character*8 :: passfilename
    integer i,j,k,l,m,n

    if(.not.allocated(tempR))    allocate(tempR(Rn3,Rn2,Rn1));          tempR=0.
    if(.not.allocated(tempC))    allocate(tempC(Cn3,Cn2,Cn1));          tempC=0.
    if(.not.allocated(uNM1))     allocate(uNM1(Rn3,Rn2,Rn1,nAtom));     uNM1=(0.,0.)
    if(.not.allocated(uNM2))     allocate(uNM2(Rn3,Rn2,Rn1,nAtom));     uNM2=(0.,0.)
    if(.not.allocated(uNM3))     allocate(uNM3(Rn3,Rn2,Rn1,nAtom));     uNM3=(0.,0.)
    if(.not.allocated(q0Gr))     allocate(q0Gr(Rn3,Rn2,Rn1));           q0Gr=0.
    if(.not.allocated(fExpNM))   allocate(fExpNM(Rn3,Rn2,Rn1,nAtom));   fExpNM=(0.,0.)
    if(.not.allocated(fExpN))    allocate(fExpN(Rn3,Rn2,Rn1,nAtom));    fExpN=(0.,0.)
    if(.not.allocated(fU1ExpN))  allocate(fU1ExpN(Rn3,Rn2,Rn1,nAtom));  fU1ExpN=(0.,0.)
    if(.not.allocated(fU2ExpN))  allocate(fU2ExpN(Rn3,Rn2,Rn1,nAtom));  fU2ExpN=(0.,0.)
    if(.not.allocated(fU3ExpN))  allocate(fU3ExpN(Rn3,Rn2,Rn1,nAtom));  fU3ExpN=(0.,0.)
    if(.not.allocated(fExpNRk))  allocate(fExpNRk(Cn3,Cn2,Cn1,nAtom));  fExpNRk=(0.,0.)
    if(.not.allocated(fU1ExpNRk))allocate(fU1ExpNRk(Cn3,Cn2,Cn1,nAtom));fU1ExpNRk=(0.,0.)
    if(.not.allocated(fU2ExpNRk))allocate(fU2ExpNRk(Cn3,Cn2,Cn1,nAtom));fU2ExpNRk=(0.,0.)
    if(.not.allocated(fU3ExpNRk))allocate(fU3ExpNRk(Cn3,Cn2,Cn1,nAtom));fU3ExpNRk=(0.,0.)
    if(.not.allocated(fExpNIk))  allocate(fExpNIk(Cn3,Cn2,Cn1,nAtom));  fExpNIk=(0.,0.)
    if(.not.allocated(fU1ExpNIk))allocate(fU1ExpNIk(Cn3,Cn2,Cn1,nAtom));fU1ExpNIk=(0.,0.)
    if(.not.allocated(fU2ExpNIk))allocate(fU2ExpNIk(Cn3,Cn2,Cn1,nAtom));fU2ExpNIk=(0.,0.)
    if(.not.allocated(fU3ExpNIk))allocate(fU3ExpNIk(Cn3,Cn2,Cn1,nAtom));fU3ExpNIk=(0.,0.)
    if(.not.allocated(AmpA))     allocate(AmpA(Cn3,Cn2,Cn1));           AmpA=(0.,0.)
    if(.not.allocated(AmpB))     allocate(AmpB(Cn3,Cn2,Cn1));           AmpB=(0.,0.)
    if(.not.allocated(IA))       allocate(IA(Cn3,Cn2,Cn1));             IA=0.
    if(.not.allocated(IB))       allocate(IB(Cn3,Cn2,Cn1));             IB=0.

    call MPI_Bcast(eAvg,3,MPI_real8,0,MPI_Comm_World,ierr)

    G1 = G10 - qC10*eAvg(1) - qC20*eAvg(6) - qC30*eAvg(5)
    G2 = G20 - qC10*eAvg(6) - qC20*eAvg(2) - qC30*eAvg(4)
    G3 = G30 - qC10*eAvg(5) - qC20*eAvg(4) - qC30*eAvg(3)
!test    if(rank==0) write(*,*) "eAvg =",eAvg(1),eAvg(2),eAvg(3)

    qC1 = qC10 + G1 - G10
    qC2 = qC20 + G2 - G20
    qC3 = qC30 + G3 - G30

    mqk1_out = dmqk1 + qC1
    mqk2_out = dmqk2 + qC2
    mqk3_out = dmqk3 + qC3

    q0Gr = 0.
    do i=1,Rn1;  q0Gr(:,:,i)=q0Gr(:,:,i)+(qC1-G1)*(i+lstartR)*dx;  enddo
    do i=1,Rn2;  q0Gr(:,i,:)=q0Gr(:,i,:)+(qC2-G2)*i*dy;  enddo
    do i=1,Rn3;  q0Gr(i,:,:)=q0Gr(i,:,:)+(qC3-G3)*i*dz;  enddo

    fExpN = 0.d0
    fU1ExpN = 0.d0
    fU2ExpN = 0.d0
    fU3ExpN = 0.d0

    do n = 1,nAtom
      do m = 1,nPhase
        uNM1(:,:,:,n) = u1(:,:,:)
        uNM2(:,:,:,n) = u2(:,:,:)
        uNM3(:,:,:,n) = u3(:,:,:)
        do k = 1,nStruc
          uNM1(:,:,:,n) = uNM1(:,:,:,n) + su1(k,m,n)*os(k,:,:,:)
          uNM2(:,:,:,n) = uNM2(:,:,:,n) + su2(k,m,n)*os(k,:,:,:)
          uNM3(:,:,:,n) = uNM3(:,:,:,n) + su3(k,m,n)*os(k,:,:,:)
        enddo
        fExpNM(:,:,:,n) = cdexp( -(0,1)* (qC1*uNM1(:,:,:,n) + qC2*uNM2(:,:,:,n) + qC3*uNM3(:,:,:,n) + q0Gr) ) *fAtom(m,n)*op(m,:,:,:)*region
          fExpN(:,:,:,n) =   fExpN(:,:,:,n) + fExpNM(:,:,:,n)
        fU1ExpN(:,:,:,n) = fU1ExpN(:,:,:,n) + uNM1(:,:,:,n) * fExpNM(:,:,:,n)
        fU2ExpN(:,:,:,n) = fU2ExpN(:,:,:,n) + uNM2(:,:,:,n) * fExpNM(:,:,:,n)
        fU3ExpN(:,:,:,n) = fU3ExpN(:,:,:,n) + uNM3(:,:,:,n) * fExpNM(:,:,:,n)
      enddo
    enddo

    do n = 1,nAtom
      tempR = real(  fExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);    fExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(  fExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);    fExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = real(fU1ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU1ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(fU1ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU1ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = real(fU2ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU2ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(fU2ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU2ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = real(fU3ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU3ExpNRk(:,:,:,n) = tempC/(nx*ny*nz)
      tempR = imag(fU3ExpN(:,:,:,n));  call mupro_fft_forward(tempR,tempC);  fU3ExpNIk(:,:,:,n) = tempC/(nx*ny*nz)
    enddo

    AmpA = (0.,0.)
    AmpB = (0.,0.)
    do n=1,nAtom
      AmpA = AmpA + fA(:,:,:,n) * (                          fExpNRk(:,:,:,n)  + (0,1)*          fExpNIk(:,:,:,n)   &
                                    -(0,1)* (mqkA1*       fU1ExpNRk(:,:,:,n)  + mqkA2*       fU2ExpNRk(:,:,:,n)  + mqkA3*       fU3ExpNRk(:,:,:,n)) &
                                    +       (mqkA1*       fU1ExpNIk(:,:,:,n)  + mqkA2*       fU2ExpNIk(:,:,:,n)  + mqkA3*       fU3ExpNIk(:,:,:,n))  )
      AmpB = AmpB + fB(:,:,:,n) * (                 dconjg(  fExpNRk(:,:,:,n)) + (0,1)* dconjg(  fExpNIk(:,:,:,n))  &
                                    -(0,1)* (mqkB1*dconjg(fU1ExpNRk(:,:,:,n)) + mqkB2*dconjg(fU2ExpNRk(:,:,:,n)) + mqkB3*dconjg(fU3ExpNRk(:,:,:,n)))&
                                    +       (mqkB1*dconjg(fU1ExpNIk(:,:,:,n)) + mqkB2*dconjg(fU2ExpNIk(:,:,:,n)) + mqkB3*dconjg(fU3ExpNIk(:,:,:,n))) )
    enddo

    AmpA = AmpA * kfExpA
    AmpB = AmpB * kfExpB

    IA = max(cdabs(AmpA)**2,1.d-50)
    IB = max(cdabs(AmpB)**2,1.d-50)

    call ArrayFourierToRegular(IA,IB,I_out,trans)

    deallocate(tempR,tempC)
    deallocate(q0Gr)
    deallocate(uNM1,uNM2,uNM3,fExpNM)
    deallocate(fExpN,fU1ExpN,fU2ExpN,fU3ExpN)
    deallocate(fExpNRk,fU1ExpNRk,fU2ExpNRk,fU3ExpNRk)
    deallocate(fExpNIk,fU1ExpNIk,fU2ExpNIk,fU3ExpNIk)
    deallocate(AmpA,AmpB,IA,IB)

  end subroutine



  subroutine ArrayFourierToRegular(ArrayA_in,ArrayB_in,Array_out,trans_in)

    use mod_interface_diffraction
    use Size_FFT_Para

    implicit none

    real*8,intent(IN),dimension(:,:,:) :: ArrayA_in,ArrayB_in
    real*8,intent(OUT),dimension(:,:,:) :: Array_out
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
!test    write(*,567),rank,Cn1Sum,Cn1,Cn2,Cn3,lstartCAll(rank),lstart3

    if(.not.allocated(ArrayA_All))    allocate(ArrayA_All(Cn3,Cn2,Cn1Sum));  ArrayA_All=0.
    if(.not.allocated(ArrayB_All))    allocate(ArrayB_All(Cn3,Cn2,Cn1Sum));  ArrayB_All=0.
    if(.not.allocated(Array_All))     allocate(Array_All(nz,ny,nx));         Array_All=0.
    if(.not.allocated(Array_out_All)) allocate(Array_out_All(nz,ny,nx));     Array_out_All=0.

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

!test    write(*,*) "Array_All size =", size(Array_All,1), size(Array_All,2), size(Array_All,3)
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
