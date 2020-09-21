#============================================================================================
#	Modo de empleo:  make [opcion]
#
#       Opciones: 
#
#                 => Compila el codigo sph y genera el ejecutable sph.e
#         clean   => Limpia los restos de las compilaciones y las ejecuciones
#=============================================================================================


#=============================
# Code structure
#=============================

AQUI    = work
GRAFS   = work/graficos
#RES     = $(DIR)/Outdata
#UTIL    = $(DIR)/Utilities
#GRAFS   = $(DIR)/Grafs
#IN      = $(DIR)/Indata
#OUT     = $(DIR)/Outdata

#=============================
# Main code and executable 
#=============================

SOURCE  = main.f90      
TARGET  = sph.e
#INCL    =-I/usr/local/include 
#LIBS    = $(INCL) -L/usr/local/lib  -L/opt/intel/composer_xe_2013_sp1.0.080/compiler/lib/intel64/ 
LIBS    = 

#=============================
# Define objects
#=============================

NUCL	= jjose2.o burn.o
#FOBJS   = mod_essentials.o mod_parameters.o mod_functions.o mod_commons.o             \#
FOBJS   = mod_parameters.o mod_functions.o mod_commons.o  mod_EOS.o\
          startout.o inparams.o indata.o                                              \
	  iter_rhoh.o treeconst.o indexx.o indexxi2.o                                 \
	  EOS.o EOS0.o helmholtz.o degenerate.o degenerate2.o $(NUCL)                 \
          relax.o relax_frame.o approach.o separate.o                                 \
          diagnostics.o outdata.o outdata_ascii.o outchem.o gettime.o energy.o                        \
	  kickstart.o predcorr.o noninertial.o forces.o varydt.o hydro_rs.o           \
          predictor.o corrector.o sort.o\
	  layer.o norm_layer.o putlayer.o                 \
          genera_single.o genera_bin.o sph.o analysis.o plot.o main.o
OBJS    = $(FOBJS)

#=============================
#Precompiler options
#=============================
ifeq ($(debug), yes)
    FPPFLAGS += -Ddebug -fpe0 -check all -debug all -traceback
#    FPPFLAGS += -Ddebug -traceback
endif
ifeq ($(openmp), yes)
   PARFLAG   = -openmp
   FPPFLAGS += -Dopenmp
endif
ifeq ($(Helium),yes)
   FPPFLAGS += -DHelium
endif
#PARFLAG   = -openmp
#FPPFLAGS += -Dopenmp
#FPPFLAGS += -Dopenmp -DHelium

#=============================
#Define Compilers
#=============================

ifeq ($(MPI), yes)
FC	= mpifort -O3 -r8 -ipo -mcmodel medium -shared-intel -DMPI $(PARFLAG) $(FPPFLAGS)
else
FC      = gfortran $(PARFLAG) $(FPPFLAGS)
endif

#=============================
#Linking code
#=============================

$(TARGET):  $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS)


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
