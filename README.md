
# Satellite Simulator and Sandbox for Cloud Observation and Modelling


**A satellite simulator and retrieval sandbox tool for cloud studies.**

S3COM is a tool at the interface of high-resolution modeling and satellite observations. It aims to make cloud studies a little easier by
- Simulating realistic satellite observations and cloud products from model outputs
- Quantifying the sensitivity of radiative quantities to cloud parameters (in progress)
- Assisting the development of retrieval algorithms using output fields from high-resolution models (in progress)

## Install

S3COM can be installed via its makefile:
```bash
git clone https://github.com/odrans/S3COM.git
cd S3COM
make
```

<details>
  <summary> <b>Click here</b> for more details on S3COM dependencies and environment.</summary>

**Dependencies**

The following dependencies are needed and require path adjustments in the Makefile:
- [NetCDF4](https://www.unidata.ucar.edu/software/netcdf/) (C and Fortran)
  - **Makefile path:** `PATH_NETCDF_C`, `PATH_NETCDF_F`
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
  - **Makefile path:** `PATH_HDF5`
- [RTTOV](https://nwp-saf.eumetsat.int/site/software/rttov)
  - RTTOV v13 is the main radiative transfer code in S3COM.
  - **Makefile path:** `PATH_RTTOV`
  - **Notes:**
    - RTTOV can be installed following these [instructions](https://nwp-saf.eumetsat.int/site/software/rttov/rttov-v13/#Installing_RTTOV_v132). Note that it requires NetCDF4 and HDF5 to be explicitely linked in its `build/Makefile.local` file before running `rttov_compile.sh`.
    - Depending on the spectrum of simulated instruments, adequate [coefficient files](https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/coefficient-download/) and [BRDF / emissivity atlases](https://nwp-saf.eumetsat.int/site/software/rttov/download/#Emissivity_BRDF_atlas_data) should be downloaded and included in the RTTOV repository.
    - Following RTTOV recommendations, it is advised to raise the system stack size: `ulimit -s unlimited`.

All Makefile paths refer to base repositories giving access to libraries, objects and binaries.

**Environment**
  
`PATH_S3COM` should be set in Makefile to indicate where S3COM is installed.

[The environment section](docs/environment.md) provides advised settings on specific supercomputers.

</details>

Note: Run `make clean` to remove all built binaries, libraries and objects.

## Usage

S3COM requires the path to a namelist file as main argument:

```bash
./s3com config_default.nml
```
The namelist contains information on input and output files as well as important options to run S3COM. Refer to the [namelist section](docs/namelist.md) for a detailed description.


**S3COM outputs** are 3 files containing satellite simulations, retrievals and atmospheric data. See the [output section](docs/outputs.md).


## Current limitations

- The current version is only configured to handle simulations from the **ICON** model. Some adjustments will be necessary to use outputs from another model. 
- S3COM only simulates measurements from passive remote-sensing sensors.
- Polarization and 3D effects are not included.

## R package

The [Rs3com](https://github.com/odrans/Rs3com) package was developed to conveniently create namelist files, run the algorithm, read its input and output and create basic figures.

## Caution

S3COM does not account for sub-grid variability of cloud properties. It should be used on atmospheric model simulations with a spatial resolutions similar or higher than that of the selected satellite instruments, ideally from CRMs or LES models. It is advised against using S3COM on GCM outputs; for such use we advise more dedicated satellite simulators, such as [COSPv2](https://github.com/CFMIP/COSPv2.0). 

## License

S3COM is available under a BSD 3-clause license. See [LICENSE](https://github.com/odrans/S3COM/blob/main/LICENSE) for details.

