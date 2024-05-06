
include ./Makefile.config

# level-1 modules: self-contained
# --level-2 modules: use level-1 modules and are used by level-3 modules
#   --level-3 modules: use both level-1 and level-2 modules
#
OBJS = $(DIR_OBJ)/global.o    \
       $(DIR_OBJ)/parser.o    \
       $(DIR_OBJ)/datatype.o    \
       $(DIR_OBJ)/mmi_mpi.o   \
       $(DIR_OBJ)/modeling.o  \
       $(DIR_OBJ)/record_wave_boundary_uw.o  \
       $(DIR_OBJ)/string.o    \
       $(DIR_OBJ)/su.o        \
       $(DIR_OBJ)/math.o      \
       $(DIR_OBJ)/crosco.o      \
       $(DIR_OBJ)/free_surface.o      \
       $(DIR_OBJ)/io.o        \
       $(DIR_OBJ)/smooth_slow.o  \
       $(DIR_OBJ)/pml.o       \
       $(DIR_OBJ)/source.o    \
       $(DIR_OBJ)/xcorr_calc_functions.o \
       $(DIR_OBJ)/gradient_precon.o \
       $(DIR_OBJ)/staggered_elastic_pml_time_step.o \
       $(DIR_OBJ)/staggered_elastic_kernel.o \
       $(DIR_OBJ)/modeling_elastic_main.o \
       $(DIR_OBJ)/fwi_elastic_subroutine.o \
       $(DIR_OBJ)/fwi_elastic_main.o \

all: dir main_fi_elas

dir:
	@-test -d $(DIR_MOD) || mkdir $(DIR_MOD)
	@-test -d $(DIR_OBJ) || mkdir $(DIR_OBJ)
	@-test -d $(DIR_EXE) || mkdir $(DIR_EXE)

main_fi_elas: $(OBJS) $(DIR_OBJ)/main_fi_elas.o
	$(FCC) -o $@ $^ $(LDFLAGS)
	mv $@ $(DIR_EXE)

bin2su.x: $(OBJS) $(DIR_OBJ)/bin2su.o
	$(FCC) -o $@ $^ $(LDFLAGS)
	mv $@ $(DIR_EXE)

$(DIR_OBJ)/%.o: %.c
	$(CC) -o $@ $< $(CCFLAGS)

$(DIR_OBJ)/%.o: %.f90
	$(FCC) -o $@ $< $(F90FLAGS)

clean:
	rm -f $(DIR_EXE)/* $(DIR_OBJ)/*.o* $(DIR_MOD)/*.mod;



       ##$(DIR_OBJ)/m_fft_utils.o    \
       ##$(DIR_OBJ)/fft.o    \

