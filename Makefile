#============================================================================================
#	Modo de empleo:  make [opcion]
#
#       Opciones: 
#
#                 => Compila el codigo sph y genera el ejecutable sph.e
#         clean   => Limpia los restos de las compilaciones y las ejecuciones
#=============================================================================================

#=============================
# Main code and executable 
#=============================

SOURCE  = main.f90      
TARGET  = sph.e

#=============================
# Define objects
#=============================

#NUCL	= jjose2.o burn.o
OBJS  = mod_parameters.o mod_functions.o mod_commons.o  mod_eos.o                  \
        startout.o inparams.o indata.o iter_rhoh.o treeconst.o indexx.o indexxi2.o \
        eos.o eos0.o helmholtz.o degenerate.o degenerate2.o $(NUCL)                \
        relax.o relax_frame.o approach.o separate.o                                \
        diagnostics.o outdata.o outchem.o gettime.o energy.o                       \
	      kickstart.o predcorr.o noninertial.o forces.o varydt.o hydro_rs.o          \
        predictor.o corrector.o sort.o layer.o norm_layer.o putlayer.o             \
        genera_single.o genera_bin.o sph.o analysis.o plot.o main.o

#=============================
#Precompiler options
#=============================
ifeq ($(debug), yes)
   FPPFLAGS += -Ddebug -fpe0 -check all -debug all -traceback
endif
ifeq ($(openmp), yes)
   PARFLAG   = -qopenmp
   FPPFLAGS += -Dopenmp
endif
ifeq ($(Helium),yes)
   FPPFLAGS += -DHelium
endif
ifeq ($(global), yes)
   FPPFLAGS  += -Dglobal
endif
ifeq ($(oldformat),yes)
   FPPFLAGS += -Doldformat
endif

#=============================
#Default options
#=============================

#PARFLAG   = -qopenmp
#FPPFLAGS += -Dopenmp

#=============================
#Define Compilers
#=============================

ifeq ($(MPI), yes)
FC	= mpiifort -O3 -xCORE-AVX512 -mcmodel=medium -warn uninitialized -warn truncated_source\
			-warn interfaces -nogen-interfaces -DMPI $(PARFLAG) $(FPPFLAGS)
else
#FC	= ifort -O3 -r8 -xCORE-AVX512 -mcmodel=medium -warn uninitialized -warn truncated_source\
			-warn interfaces -nogen-interfaces $(PARFLAG) $(FPPFLAGS)
FC	= gfortran -Dopenmp -fopenmp -fdefault-real-8
endif


#=============================
#Linking code
#=============================

$(TARGET):  $(OBJS)
#	$(FC) -o $@ $(OBJS) $(LIBS)
	$(FC) -o sph.e $(OBJS) -Dopenmp -fopenmp -fdefault-real-8

#=============================
#Compile code
#=============================

# Fortran 90 subroutines
%.o : %.f90
	$(FC) -c  $<

%.o : %.F90
	$(FC) -c  $<

%.o : %.f
	$(FC) -c  $<

#================================
#Clean compilation 
#================================

clean:	
	rm *.o *.mod 

#================================
#Start simulation
#================================

run:
	nohup mpirun -np 4 --hostfile mylist Gabri.e &

#================================
#Clean ALL simulation results
#================================

new:
	rm -f fort* nohup.out
	rm -f *.out degdt h*0

#================================
#Plot results
#================================

plot: 
	cp $(File) aux
	cp $(CFile) caux 
	ifort $(SUBS)/dibujos_mm.f90 -o plots
	./plots
	ifort $(SUBS)/out_energy.f -o energs
	./energs
	rm plots energs
	mv aux caux $(GRAFS)
	plotmtv $(GRAFS)/en_picture $(GRAFS)/perf_* 
