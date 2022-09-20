Satellite Simulator and Sandbox for Cloud Observation and Modelling (S3COM)
================================================================================

S3COM is a satellite simulator and retrieval sandbox tool for cloud studies. It provides satellite observations and cloud products based on high-resolution modelling outputs.

The main goals of S3COM are:
- providing realistic satellite measurements consistent with model outputs
- computing the sensitivity of radiative quantities to cloud parameters
- evaluating and improving retrieval algorithms using output fields from high-resolution models

How to use
----------

- Section in progress

Environment on DKRZ/Levante: module purge && module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0 openmpi/4.1.2-intel-2021.5.0 netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0 netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0 hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0

Edit the pathnames to your input ("fname_in") and output ("fname_out") netcdf files in config.nml

Compiling & Requirements
------------------------

[RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov) v13.1 is the main radiative transfer algorithm used in S3COM. It needs to be installed on your system and its libraries linked in the Makefile via the `RTTOV_PATH` variable. Please refer to the RTTOV documentation. S3COM has not been tested for other versions of RTTOV. Following RTTOV recommendations, it is advised to raise the system stack size, e.g.: `ulimit -s unlimited` 

NetCDF4 (C and Fortran) and HDF5 libraries are required. It is advised to link them by editing `PATH_NCDF_C_LIB`, `PATH_NCDF_LIB`, `PATH_NCDF_INC` and `PATH_HDF5_LIB` accordingly in the Makefile. 

The `basedir` variable must also be edited in the Makefile to indicate the repository where S3COM is installed.

`make clean` cleans the S3COM repositories from binaries and installed libraries and modules. `make install` (or `make`) compiles the code to create the `s3com` binary.

Current limitations
-------------------

- The current version is only configured to handle simulations from the ICON model. Some adjustments will be necessary to use outputs from another model. 
- S3COM only simulates measurements from passive remote-sensing sensors.
- Polarization and 3D effects are not yet included.

Caution
-------

S3COM does not account for sub-grid variability of cloud properties. It should be used on atmospheric model simulations with a spatial resolutions similar or higher than that of the selected satellite instruments, ideally from CRMs or LES models. It is advised against using S3COM on GCM outputs; for such use we advise more dedicated satellite simulators, such as [COSPv2](https://github.com/CFMIP/COSPv2.0). 

License
------
S3COM is available under a BSD 3-clause license.
Please see [LICENSE](LICENSE) for details.
