# How to install S3COM

This guide will help you install the required dependencies for S3COM: HDF5, NetCDF4, and RTTOV. 

## RTTOV

Follow the instructions on the [EUMETSAT webiste](https://nwp-saf.eumetsat.int/site/software/rttov/rttov-v13) on how to install RTTOV.

Briefly:
- Create an account if you haven't already
- Register for RTTOV v13 (version currently supported for S3COM)
- download rttov{version}.tar.xz and untar it, e.g.
```bash
mkdir rttov132
cd rttov132
tar xf rttov132.tar.xz
```
- download the gas, aerosol and cloud coefficient files:
```bash
cd rtcoef_rttov13
./rttov_coef_download.sh
```
  - download all coefficients if unsure of what you need.
- download the necessary [emissivity atlases](https://nwp-saf.eumetsat.int/site/software/rttov/download)
  - to be stored in `emis_data`
  - Example of script for downloading the entire emissivity atlas dataset:
      ```bash
      #!/bin/bash

      urls=(
          "https://nwp-saf.eumetsat.int/downloads/emis_data/uw_ir_emis_atlas_hdf5.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/camel_emis_atlas_v01r02.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/uw_ir_emis_atlas_angcorr_hdf5.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/camel_clim_emis_atlas_pchsr.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/camel_clim_emis_atlas_jan-mar.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/camel_clim_emis_atlas_apr-jun.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/camel_clim_emis_atlas_jul-sep.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/camel_clim_emis_atlas_oct-dec.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/uw_ir_emis_atlas_angcorr_hdf5.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/telsem2_mw_atlas.tar.bz2"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/cnrm_mwemis_amsu_mhs_data.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/cnrm_mwemis_atms_data.tar"
          "https://nwp-saf.eumetsat.int/downloads/emis_data/cnrm_mwemis_ssmis_data.tar"
      )

      for url in "${urls[@]}"; do
          filename=$(basename "$url")
          wget "$url" && tar xf "$filename" && rm "$filename"
      done
      ```
- Simiarly download the BRDF atlas
  - to be stored in `brdf_data`
  - As of v13.2, a single file to download:
    ```bash 
       wget https://nwp-saf.eumetsat.int/downloads/brdf_data/cms_brdf_atlas_hdf5.tar
       tar xf cms_brdf_atlas_hdf5.tar
    ```
- edit the following variables in `build/Makefile.local`:
  - `HDF5_PREFIX`: path to `hdf5`. Also uncomment `FFLAGS_HDF5` and `LDFLAGS_HDF5`. Make sure that the hdf5hl library is `-lhdf5hl_fortran`, and **not** `-lhdf5_hl_fortran` (bug in RTTOV v13).
  - `NETCDF_PREFIX`: path to `netcdf-fortran` (also uncomment `FFLAGS_NETCDF` and `LDFLAGS_NETCDF`)
  - We advise to use LAPACK libraries contained in RTTOV (i.e. keep `LAPACK_PREFIX` commented)
- Compile RTTOV using the interactive script:
```bash
cd src
../build/rttov_compile.sh
```
  - S3COM has been tested with ifort and gfortran. Use `gfortran-openmp` and `ifort-openmp` for OpenMP support and multithread use of S3COM.






