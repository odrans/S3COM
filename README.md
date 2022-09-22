
# Satellite Simulator and Sandbox for Cloud Observation and Modelling (S3COM)

**A satellite simulator and retrieval sandbox tool for cloud studies.**

**S3COM** aims to make cloud studies a little easier by
- Providing realistic satellite measurements and cloud products consistent with model outputs
- Computing the sensitivity of radiative quantities to cloud parameters
- Assisting the development of retrieval algorithms using output fields from high-resolution models

## How to use

Important parameters to run **S3COM** are stored in a namelist file. The default namelist is `config_default.nml`. Some path must be edited before running S3COM for the first time. See the [namelist section](namelist.md) for details.

## Environment & Compiling

[RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov) v13.1 is the main radiative transfer algorithm used in **S3COM**. It needs to be installed on your system and its libraries linked in the Makefile via the `RTTOV_PATH` variable. Please refer to the RTTOV documentation. **S3COM** has not been tested for other versions of **RTTOV**. Following RTTOV recommendations, it is advised to raise the system stack size, e.g.: `ulimit -s unlimited` 

**NetCDF4** (C and Fortran) and **HDF5** libraries are required. It is advised to link them by editing `PATH_NCDF_C_LIB`, `PATH_NCDF_LIB`, `PATH_NCDF_INC` and `PATH_HDF5_LIB` accordingly in the Makefile. 

The `basedir` variable must also be edited in the Makefile to indicate the repository where **S3COM** is installed.

 `make install` (or just `make`) compiles the code to create the `s3com` binary. `make clean` cleans all repositories.

Advised environments for specific supercomputers can be found [here](Environment.md).

## Current limitations

- The current version is only configured to handle simulations from the **ICON** model. Some adjustments will be necessary to use outputs from another model. 
- **S3COM** only simulates measurements from passive remote-sensing sensors.
- Polarization and 3D effects are not included.

## R package

The [Rs3com](https://github.com/odrans/Rs3com) package was developed to conveniently create namelist files, run the algorithm, read its input and output and create basic figures.

## Caution

**S3COM** does not account for sub-grid variability of cloud properties. It should be used on atmospheric model simulations with a spatial resolutions similar or higher than that of the selected satellite instruments, ideally from CRMs or LES models. It is advised against using **S3COM** on GCM outputs; for such use we advise more dedicated satellite simulators, such as [COSPv2](https://github.com/CFMIP/COSPv2.0). 

## License

**S3COM** is available under a BSD 3-clause license.
Please see [LICENSE](LICENSE) for details.
