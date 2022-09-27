
# Satellite Simulator and Sandbox for Cloud Observation and Modelling | S3COM

**A satellite simulator and retrieval sandbox tool for cloud studies.**

S3COM is a tool at the interface of high-resolution modeling and satellite observations. It aims to make cloud studies a little easier by
- Computing realistic satellite observations and cloud products from model outputs
- Quantifying the sensitivity of radiative quantities to cloud parameters (in progress)
- Assisting the development of retrieval algorithms using output fields from high-resolution models (in progress)

## Install

S3COM can be installed via its makefile:
```bash
git clone git@github.com:odrans/S3COM.git
cd S3COM
make
```

<details>
  <summary> <b>Click here</b> for more details on S3COM dependencies and environment.</summary>

### Dependencies

The following dependencies are required and should be adjusted in the Makefile:
- [RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov)
  - RTTOV v13.1 is the main radiative transfer code in S3COM. Please carefuly refer its documentation for the installation.
  - **Makefile path:** `PATH_RTTOV`
  - **Notes:** 
    - Following RTTOV recommendations, it is advised to raise the system stack size: `ulimit -s unlimited`.
    - Depending on the spectrum of simulated instruments, adequate [coefficient files](https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/coefficient-download/) and [BRDF / emissivity atlases](https://nwp-saf.eumetsat.int/site/software/rttov/download/#Emissivity_BRDF_atlas_data) should be downloaded and added to the RTTOV repository.
- [NetCDF4](https://www.unidata.ucar.edu/software/netcdf/) (C and Fortran)
  - **Makefile path:** `PATH_NETCDF_C`, `PATH_NETCDF_F`
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
  - **Makefile path:** `PATH_HDF5`
  
All Makefile paths refer to base repositories, giving access to libraries, objects and binaries.

### Environment  
  
`PATH_S3COM` should be set in Makefile to indicate where S3COM is installed.

[The environment section](Environment.md) provides advised settings on specific supercomputers.

</details>

Note: Run `make clean` to remove all built binaries, libraries and objects.

## Usage

S3COM requires the path to a namelist file as main argument:

```bash
./s3com config_default.nml
```
The namelist contains information on input and output files as well as important options to run S3COM. Refer to the [namelist section](namelist.md) for a detailed description.


**S3COM outputs** are 3 files containing satellite simulations, retrievals and atmospheric data. See the [output section](output.md).


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
