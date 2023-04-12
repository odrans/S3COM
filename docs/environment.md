# Environments and job scripts for specific servers

S3COM can be executed localy as well as on supercomputers. Note that all libraries, including RTTOV, must be compiled with openMP or similar to allow for parallelization (see `rttov_nthreads` option in namelist).

## Levante

Example of modules to load to run RTTOV with multiple threads. 

With intel compilers:
```bash
module purge
module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0 \
            openmpi/4.1.2-intel-2021.5.0 \
            netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0 \
            netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0 \
            hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0 \
            eccodes/2.21.0-intel-2021.5.0
```

With gcc compilers:
```bash
module load gcc/11.2.0-gcc-11.2.0 \
            openmpi/4.1.2-gcc-11.2.0 \
            parallel-netcdf/1.12.2-openmpi-4.1.2-gcc-11.2.0 \
            netcdf-fortran/4.5.3-openmpi-4.1.2-gcc-11.2.0 \
            netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0 \
            hdf5/1.12.1-openmpi-4.1.2-gcc-11.2.0 \
            eccodes/2.21.0-gcc-11.2.0
```

Job script (edit project number):

```bash
#!/bin/bash
#SBATCH -J S3COM                                  # Specify job name
#SBATCH -p shared                                 # Use partition shared
#SBATCH --nodes=1                                 # Number of requested nodes
#SBATCH --ntasks-per-node=48                      # Tasks per nodes
#SBATCH -A PROJECT_NUMBER                         # Charge resources on this project account
#SBATCH -o s3com.o%j                              # File name for standard and error output

# Environment settings to run a MPI parallel program
module purge
module load intel-oneapi-compilers
module load openmpi/4.1.2-intel-2021.5.0
module load hdf5/1.12.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-cxx/4.2-openmpi-4.1.2-intel-2021.5.0
module load parallel-netcdf/1.12.2-openmpi-4.1.2-intel-2021.5.0
module load libaec/1.0.5-intel-2021.5.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack-levante/hdf5-1.12.1-tvymb5/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib

# Limit stacksize ... adjust to your programs need and core file size
ulimit -s unlimited
ulimit -c 0

# Execute OpenMP programs
srun ./s3com
```
