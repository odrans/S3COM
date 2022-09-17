Satellite Simulator and Sandbox for Cloud Observation and Modelling (S3COM)
================================================================================

S3COM is a satellite simulator and sandbox tool for cloud studies. It simulates satellite observations consistent with high-resolution modelling outputs and can also use these simulations to provide cloud retrievals. It currently focuses on passive satellite observations and products. 

How to use
--------------------------------------------------------------------------------

g-- Section in progress

Environment on DKRZ/Levante:
module purge && module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0 openmpi/4.1.2-intel-2021.5.0 netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0 netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0 hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0

Edit the pathnames to your input ("fname_in") and output ("fname_out") netcdf files in config.nml

In the Makefile:
- edit the paths to:
   - netcdf-c library ("path_NCDF_C_LIB")
   - netcdf-fortran library ("path_NCDF_INC" and "path_NCDF_LIB")
   - HDF5 library ("HDF5_LIB")
- edit the path to rttov v13.1 ("RTTOV_PATH").

After modifications in the code, type the following commands in the terminal being in "/home/ML_RTTOV_CDNC/":
- make clean (cleaning the code)
- make or make install (compiling the code)
- ulimit -s unlimited
- ./s3com (running the code)

License
--------------------------------------------------------------------------------
S3CM is available under a BSD 3-clause license.
Please see [LICENSE](LICENSE) for details.
