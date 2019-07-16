TTM=pot_ttm

FILES=four1 pflush fsiesta_pipes system_mod lun_management geometry_calculator diffusion_calculator dans_timer consts $(TTM)/math spectral_properties Nose_Hoover NormalModes $(TTM)/pot_mod Langevin main_stuff estimators $(TTM)/find_neigh $(TTM)/nasa_mod $(TTM)/nasa $(TTM)/ewald $(TTM)/pot_spc $(TTM)/pot_ttm $(TTM)/potential $(TTM)/smear dip_ttm  input_file_reader InputOutput  force_calc MD PIMD

OBJS=$(addsuffix .o, $(FILES))

#VPATH=/home/dan/Dropbox/RESEARCH/SIESTA_myMonomerCorrected/Src

FC=mpif90

##Note : -DFC_HAVE_FLUSH -DFC_HAVE_ABORT  flags are for SIESTA pipes support

#FFLAGS=-fpp  -O3 -C -debug -traceback
FFLAGS = -O3  -cpp -Dsiesta -Dparallel -DFC_HAVE_FLUSH

all: PIMD.x

#compilation for debugging
debug: FFLAGS += --debug --backtrace -fbounds-check
debug: PIMD.x

#compilation for profiling
profile: FFLAGS += -p
profile: PIMD.x

#serial compilation commands
serial: FC=gfortran  #ifort not working on Handy!
serial: FFLAGS = -O3  -cpp -Dsiesta -DFC_HAVE_FLUSH # --debug --backtrace -fbounds-check
serial: PIMD_serial.x


PIMD_serial.x: MPI.o $(OBJS)
	$(FC) $(FFLAGS) MPI.o  $(OBJS) -o ./PIMD_serial.x

PIMD.x: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o  ./PIMD.x

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o : %.F90
	$(FC) $(FFLAGS) -c $< -o $@

MPI.o:
	$(FC) $(FFLAGS) -c MPI.f90 -o MPI.o

clean:
	rm -rf *.o *.mod $(TTM)/*.o $(TTM)/*.mod

pristine:
	rm -rf *.o *.mod $(TTM)/*.o $(TTM)/*.mod INPUT_TMP* *.log *.xml *.x
