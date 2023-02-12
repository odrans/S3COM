#
#  S3COM
#  Copyright (c) 2022, University of Lille
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#   1. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#   3. Neither the name of the copyright holder nor the names of its contributors
#      may be used to endorse or promote products derived from this software without
#      specific prior written permission.

#  S3COM IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
#  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#  History
#  Jan 2022 - O. Sourdeval - Original version

prog = s3com

F90      = ifort
F90FLAGS = -module $(mod) -fpp -qopenmp -g -traceback -check bounds -debug extended
# F90FLAGS = -module $(mod) -fpp -qopenmp -march=core-avx2 -g -O3 -debug -traceback -check bounds
# F90FLAGS = -module $(mod) -fp-model source -qopenmp -g -O3 -debug -traceback -check bounds

# Let's start with setting up some paths
# ---------------------------------------------------------------------------------------------------------------------------------------
PATH_S3COM = $(HOME)/github/S3COM
PATH_RTTOV = /work/bb1036/rttov_share/rttov131
PATH_NETCDF_C = /sw/spack-levante/netcdf-c-4.8.1-2k3cmu
PATH_NETCDF_F = /sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g
PATH_HDF5 = /sw/spack-levante/hdf5-1.12.1-tvymb5
# ---------------------------------------------------------------------------------------------------------------------------------------

## Following up paths (do not edit!)
# ---------------------------------------------------------------------------------------------------------------------------------------
src = $(PATH_S3COM)/src
obj = $(PATH_S3COM)/obj
lib = $(PATH_S3COM)/lib
mod = $(PATH_S3COM)/mod

DIR_MAIN = $(src)/main
DIR_IO = $(src)/io
DIR_OE = $(src)/oe
DIR_UTILS = $(src)/utils
DIR_RTTOV = $(src)/rttov
DIR_MODELS = $(src)/models
DIR_CONF = $(src)/conf

PATH_NCDF_C_LIB = $(PATH_NETCDF_C)/lib
PATH_NCDF_INC = $(PATH_NETCDF_F)/include
PATH_NCDF_LIB = $(PATH_NETCDF_F)/lib

PATH_HDF5_LIB = $(PATH_HDF5)/lib
PATH_HDF5_INC = $(PATH_HDF5)/include

RTTOV_PATH       = $(PATH_RTTOV)
RTTOV_LIB_PATH   = $(RTTOV_PATH)/lib
RTTOV_INC_PATH   = $(RTTOV_PATH)/include
RTTOV_MOD_PATH   = $(RTTOV_PATH)/mod
RTTOV_LIBS       = -lrttov13_wrapper -lrttov13_mw_scatt -lrttov13_brdf_atlas -lrttov13_emis_atlas -lrttov13_other \
                   -lrttov13_parallel -lrttov13_coef_io -lrttov13_hdf -lrttov13_main -lhdf5_hl_fortran -lhdf5_hl \
                   -lhdf5_fortran -lhdf5
# -------------------------------------------------------------------------------------------------------------------------------

# List of library files that will be created
# -------------------------------------------------------------------------------------------------------------------------------
LIB_MAIN = $(lib)/libmain.a
LIB_IO = $(lib)/lib_io.a
LIB_OE = $(lib)/lib_oe.a
LIB_UTILS = $(lib)/libutils.a
LIB_CONF = $(lib)/libconf.a
LIB_RTTOVML = $(lib)/librttovml.a
LIB_MODELS = $(lib)/libmodels.a
# -------------------------------------------------------------------------------------------------------------------------------

# List of object files in each library
# -------------------------------------------------------------------------------------------------------------------------------
LIST_OBJ_CONF = $(obj)/types.o \
        $(obj)/config.o

LIST_OBJ_MAIN = $(obj)/s3com_setup.o

LIST_OBJ_MODELS = $(obj)/models_icon.o \
        $(obj)/models_nwpsaf.o \
        $(obj)/models.o

LIST_OBJ_IO = $(obj)/io_utils.o \
        $(obj)/io_verbose.o \
        $(obj)/io_namelist.o \
		$(obj)/io_nwpsaf.o \
		$(obj)/io_icon.o \
		$(obj)/io_s3com.o

# LIST_OBJ_OE = $(obj)/model_cloud.o \
#         $(obj)/oe_utils.o \
# 		$(obj)/oe_run.o

LIST_OBJ_UTILS = $(obj)/utils_fort.o \
		 $(obj)/utils_math.o \
		 $(obj)/utils_phys.o

LIST_OBJ_RTTOVML = $(obj)/rttov_utils.o \
		   $(obj)/rttov.o \
		   $(obj)/rttov_init.o \
		   $(obj)/rttov_setup.o

LIST_OBJ = $(LIST_OBJ_CONF) $(LIST_OBJ_UTILS) $(LIST_OBJ_IO) $(LIST_OBJ_RTTOVML) $(LIST_OBJ_MODELS) $(LIST_OBJ_OE) $(LIST_OBJ_MAIN)
# -------------------------------------------------------------------------------------------------------------------------------

# List of flags related to each libraries + final flag
# -------------------------------------------------------------------------------------------------------------------------------
FLAGS_NCDF = -I$(PATH_NCDF_INC) -L${PATH_NCDF_LIB} -lnetcdff -L${PATH_NCDF_C_LIB} -lnetcdf -Wl,-rpath,${PATH_NCDF_LIB} -Wl,-rpath,${PATH_NCDF_C_LIB}
FLAGS_RTTOV = -I${RTTOV_INC_PATH} -L${RTTOV_LIB_PATH} $(RTTOV_LIBS)
FLAG_HDF5= -L${PATH_HDF5_LIB} -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm -Wl,-rpath,${PATH_HDF5_LIB} -I${PATH_HDF5_INC}
##FLAGS_LOCAL = -L$(lib) -lmodels -l_io -l_oe -lrttovml -lmain -lutils
FLAGS_LOCAL = -L$(lib) -lmodels -l_io -lrttovml -lmain -lutils -lconf

FLAGS_ALL = $(FLAGS_LOCAL) $(FLAGS_RTTOV) $(FLAG_HDF5) $(FLAGS_NCDF)
# -------------------------------------------------------------------------------------------------------------------------------

# Make commands
# -------------------------------------------------------------------------------------------------------------------------------
install: $(LIST_OBJ)
	ar r $(LIB_CONF) $(LIST_OBJ_CONF)
	ar r $(LIB_UTILS) $(LIST_OBJ_UTILS)
	ar r $(LIB_MAIN) $(LIST_OBJ_MAIN)
	ar r $(LIB_OE) $(LIST_OBJ_OE)
	ar r $(LIB_MODELS) $(LIST_OBJ_MODELS)
	ar r $(LIB_IO) $(LIST_OBJ_IO)
	ar r $(LIB_RTTOVML) $(LIST_OBJ_RTTOVML)
	$(F90) $(F90FLAGS) $(DIR_MAIN)/$(prog).f90 -o $(prog) $(FLAGS_ALL)

clean:
	rm -f $(obj)/*.o $(mod)/*.mod $(lib)/*.a s3com
	rm -rf $(obj) $(mod) $(lib)
# -------------------------------------------------------------------------------------------------------------------------------

# Use prerequisite to make sure repositories exist
# -------------------------------------------------------------------------------------------------------------------------------
$(LIST_OBJ): | $(obj) $(lib) $(mod)

$(obj):
	mkdir -p $(obj)

$(lib):
	mkdir -p $(lib)

$(mod):
	mkdir -p $(mod)
# -------------------------------------------------------------------------------------------------------------------------------


## Objects for subroutines in ./src/conf
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/types.o : $(DIR_CONF)/types.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/config.o : $(DIR_CONF)/config.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------


## Objects for subroutines in ./src/s3com
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/s3com_setup.o : $(DIR_MAIN)/s3com_setup.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------

## Objects for subroutines in ./src/models
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/models_icon.o : $(DIR_MODELS)/models_icon.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/models_nwpsaf.o : $(DIR_MODELS)/models_nwpsaf.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/models.o : $(DIR_MODELS)/models.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------


## Objects for subroutines in ./src/oe
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/oe_run.o : $(DIR_OE)/oe_run.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/model_cloud.o : $(DIR_OE)/model_cloud.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/oe_utils.o : $(DIR_OE)/oe_utils.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------


## Objects for subroutines in ./src/io
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/io_s3com.o : $(DIR_IO)/io_s3com.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/io_icon.o : $(DIR_IO)/io_icon.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/io_nwpsaf.o : $(DIR_IO)/io_nwpsaf.f90
	$(F90) $(F90FLAGS) $(FLAG_HDF5) -I$(PATH_NCDF_INC) -c $< -o $@

$(obj)/io_verbose.o : $(DIR_IO)/io_verbose.f90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

$(obj)/io_namelist.o : $(DIR_IO)/io_namelist.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/io_utils.o : $(DIR_IO)/io_utils.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------

## Objects for subroutines in ./src/rttov
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/rttov_init.o : $(DIR_RTTOV)/rttov_init.f90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

$(obj)/rttov.o : $(DIR_RTTOV)/rttov.f90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

$(obj)/rttov_setup.o : $(DIR_RTTOV)/rttov_setup.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/rttov_utils.o : $(DIR_RTTOV)/rttov_utils.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------

## Objects for subroutines in ./src/utils
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/utils_math.o : $(DIR_UTILS)/utils_math.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/utils_fort.o : $(DIR_UTILS)/utils_fort.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/utils_phys.o : $(DIR_UTILS)/utils_phys.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------
