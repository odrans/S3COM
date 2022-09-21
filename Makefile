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

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
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
F90FLAGS = -module $(mod) -fpp -qopenmp -g -O0 -debug -traceback -check bounds
# F90FLAGS = -module $(mod) -fpp -qopenmp -march=core-avx2 -g -O3 -debug -traceback -check bounds
# F90FLAGS = -module $(mod) -fp-model source -qopenmp -g -O3 -debug -traceback -check bounds

## Let's start with some paths
# ---------------------------------------------------------------------------------------------------------------------------------------
basedir = $(HOME)/github/S3COM

src = $(basedir)/src
obj = $(basedir)/obj
lib = $(basedir)/lib
mod = $(basedir)/mod

PATH_main = $(src)/main
PATH_io = $(src)/io
PATH_utils = $(src)/utils
PATH_rttov = $(src)/rttov

PATH_NCDF_C_LIB = /sw/spack-levante/netcdf-c-4.8.1-2k3cmu/lib
PATH_NCDF_INC = /sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/include
PATH_NCDF_LIB = /sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib

PATH_HDF5_LIB = /sw/spack-levante/hdf5-1.12.1-tvymb5/lib

RTTOV_PATH       = /work/bb1036/rttov_share/rttov131
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
LIB_UTILS = $(lib)/libutils.a
LIB_RTTOVML = $(lib)/librttovml.a
# -------------------------------------------------------------------------------------------------------------------------------

# List of object files in each library
# -------------------------------------------------------------------------------------------------------------------------------
LIST_OBJ_MAIN = $(obj)/model_cloud.o \
		$(obj)/setup.o \
		$(obj)/oe_utils.o \
		$(obj)/oe_run.o

LIST_OBJ_IO =   $(obj)/regrid.o \
        $(obj)/io_namelist.o \
		$(obj)/read_icon.o \
		$(obj)/write_output.o

LIST_OBJ_UTILS = $(obj)/types.o \
		 $(obj)/config.o \
		 $(obj)/utils_math.o \
		 $(obj)/utils_fort.o

LIST_OBJ_RTTOVML = $(obj)/rttov_utils.o \
		   $(obj)/rttov_ml.o \
		   $(obj)/interface_rttov.o \
		   $(obj)/rttov_setup.o

LIST_OBJ = $(LIST_OBJ_UTILS) $(LIST_OBJ_RTTOVML) $(LIST_OBJ_IO) $(LIST_OBJ_MAIN) 
# -------------------------------------------------------------------------------------------------------------------------------

# List of flags related to each libraries + final flag
# -------------------------------------------------------------------------------------------------------------------------------
FLAGS_NCDF = -I$(PATH_NCDF_INC) -L${PATH_NCDF_LIB} -lnetcdff -L${PATH_NCDF_C_LIB} -lnetcdf -Wl,-rpath,${PATH_NCDF_LIB} -Wl,-rpath,${PATH_NCDF_C_LIB}
FLAGS_RTTOV = -I${RTTOV_INC_PATH} -L${RTTOV_LIB_PATH} $(RTTOV_LIBS)
# FLAG_HDF5= -L${PATH_HDF5_LIB} -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -Wl,-rpath,${PATH_HDF5_LIB}
FLAG_HDF5= -L${PATH_HDF5_LIB} -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm -Wl,-rpath,${PATH_HDF5_LIB}
FLAGS_LOCAL = -L$(lib) -l_io -lrttovml -lmain -lutils

FLAGS_ALL = $(FLAGS_LOCAL) $(FLAGS_RTTOV) $(FLAG_HDF5) $(FLAGS_NCDF)
# -------------------------------------------------------------------------------------------------------------------------------

# Make commands
# -------------------------------------------------------------------------------------------------------------------------------
install: $(LIST_OBJ)
	ar r $(LIB_UTILS) $(LIST_OBJ_UTILS)
	ar r $(LIB_MAIN) $(LIST_OBJ_MAIN)
	ar r $(LIB_IO) $(LIST_OBJ_IO)
	ar r $(LIB_RTTOVML) $(LIST_OBJ_RTTOVML)
	$(F90) $(F90FLAGS) $(PATH_main)/$(prog).f90 -o $(prog) $(FLAGS_ALL)

clean:	
	rm -f $(obj)/*.o $(mod)/*.mod $(lib)/*.a
# -------------------------------------------------------------------------------------------------------------------------------


## Objects for subroutines in ./src/main 
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/setup.o : $(PATH_main)/setup.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/model_cloud.o : $(PATH_main)/model_cloud.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/oe_run.o : $(PATH_main)/oe_run.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/oe_utils.o : $(PATH_main)/oe_utils.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------


## Objects for subroutines in ./src/io
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/write_output.o : $(PATH_io)/write_output.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/read_icon.o : $(PATH_io)/read_icon.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/io_namelist.o : $(PATH_io)/io_namelist.f90
	$(F90) $(F90FLAGS) -I $(PATH_NCDF_INC) -c $< -o $@

$(obj)/regrid.o : $(PATH_io)/regrid.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------

## Objects for subroutines in ./src/rttov
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/interface_rttov.o : $(PATH_rttov)/interface_rttov.f90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

$(obj)/rttov_ml.o : $(PATH_rttov)/rttov.f90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

$(obj)/rttov_setup.o : $(PATH_rttov)/rttov_setup.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/rttov_utils.o : $(PATH_rttov)/rttov_utils.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------

## Objects for subroutines in ./src/utils
# -------------------------------------------------------------------------------------------------------------------------------
$(obj)/types.o : $(PATH_utils)/types.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/config.o : $(PATH_utils)/config.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/utils_math.o : $(PATH_utils)/utils_math.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/utils_fort.o : $(PATH_utils)/utils_fort.f90 #$(PATH_utils)/config.f90
	$(F90) $(F90FLAGS) -c $< -o $@
# -------------------------------------------------------------------------------------------------------------------------------
