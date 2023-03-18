TARGET = Diffraction.exe
IFC = mpiifort

INTELSTATICFFTWROOT=/gpfs/group/lqc3/default/fftw

MKLROOT=/opt/global/RedHat-7/production/20220812/intel/compilers_and_libraries_2020.2.254/linux/mkl
MUPROROOT=/storage/work/t/tuy123/lib/MUPRO_1.0.6_Library_ElasBC_210321

INTELSTATICFFTWINC = -I${INTELSTATICFFTWROOT}/include/
INTELSTATICFFTWLIB = ${INTELSTATICFFTWROOT}/lib/libfftw3_mpi.a ${INTELSTATICFFTWROOT}/lib/libfftw3.a

#IFCFLAGS= -O0 -xSSE3 -C -warn -align rec16byte -align array64byte -align all -r8 -traceback -g
IFCFLAGS= -O3 -align rec16byte -align array64byte -align all -r8  -traceback

FC=$(IFC)
FCFLAGS=$(IFCFLAGS)
FFTWINC = $(INTELSTATICFFTWINC)
FFTWLIB = ${INTELSTATICFFTWLIB}
MKLLIB = ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a   ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -liomp5 -lm -ldl #-Wl,-rpath=XORIGIN/lib
MUPROLIB = ${MUPROROOT}/libpfMPI.a

MUPROINC = -I${MUPROROOT}
MKLINC = -I${MKLROOT}/include

SCR_F = diffraction_routines.f90 diffraction_main.f90

OBJS = $(SCR_F:.f90=.o) $(SCR_C:.c=.o)

${TARGET}: $(OBJS)  
	$(FC) $(FCFLAGS) -o $@ $(OBJS) $(MUPROLIB) $(FFTWLIB) $(MKLLIB) -lstdc++

%.o: %.f90  
	$(FC) $(FCFLAGS) $(MUPROINC) $(MKLINC) -c $< 

.PHONY:clean
clean:
	/bin/rm -f *.o *.mod *.exe *.dat fort.*
