
# Satellite Simulator and Sandbox for Cloud Observation and Modelling | S3COM

**A satellite simulator and retrieval sandbox tool for cloud studies.**

S3COM aims to make cloud studies a little easier by
- Providing realistic satellite measurements and cloud products consistent with model outputs
- Computing the sensitivity of radiative quantities to cloud parameters (in progress)
- Helping the development of retrieval algorithms using output fields from high-resolution models (in progress)

## How to use

S3COM can then be executed using a namelist file as an argument:
```
./s3com config_default.nml
```

The namelist contains information on input and output files as well as important options to run S3COM. See this [namelist description](output.md) for more detail.

## Environment & Compiling

S3COM can be compiled via its `Makefile` after adjusting a few dependencies:
- [**RTTOV**](https://nwp-saf.eumetsat.int/site/software/rttov) libraries must be linked via `RTTOV_PATH`. Please refer to the RTTOV documentation for its installation. S3COM is currently tested for RTTOV v13.1. 
- **NetCDF4** (C and Fortran) and **HDF5** libraries are required. It is advised to link them by editing `PATH_NCDF_C_LIB`, `PATH_NCDF_LIB`, `PATH_NCDF_INC` and `PATH_HDF5_LIB` accordingly. 
- `basedir`  must be set to the repository where S3COM is installed.

`make install` (or `make`) compiles the code and creates the `s3com` binary. `make clean` cleans all repositories. Note that, following RTTOV recommendations, it is advised to raise the system stack size: `ulimit -s unlimited`.

Please refer to [the environment section](Environment.md) for advised settings on specific supercomputers.

## Current limitations

- The current version is only configured to handle simulations from the **ICON** model. Some adjustments will be necessary to use outputs from another model. 
- S3COM only simulates measurements from passive remote-sensing sensors.
- Polarization and 3D effects are not included.

## R package

The [Rs3com](https://github.com/odrans/Rs3com) package was developed to conveniently create namelist files, run the algorithm, read its input and output and create basic figures.

## Caution

S3COM does not account for sub-grid variability of cloud properties. It should be used on atmospheric model simulations with a spatial resolutions similar or higher than that of the selected satellite instruments, ideally from CRMs or LES models. It is advised against using S3COM on GCM outputs; for such use we advise more dedicated satellite simulators, such as [COSPv2](https://github.com/CFMIP/COSPv2.0). 

## License

S3COM is available under a BSD 3-clause license.
Please see [LICENSE](LICENSE) for details.
