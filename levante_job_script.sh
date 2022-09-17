#!/bin/bash
#SBATCH -J ml_rttov_b381589                        # Specify job name
#SBATCH -p shared                                 # Use partition shared
##SBATCH -N 4                                       # Specify number of nodes (1 for serial applications!)
##SBATCH -n 7                                       # Specify max. number of tasks to be invoked
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
## SBATCH -t 01:00:00                               # Set a limit on the total run time
#SBATCH -A bb1143                                  # Charge resources on this project account
#SBATCH -o ml_rttov.o%j                            # File name for standard and error output
#SBATCH --mail-user=alexandre.simeon@univ-lille.fr # Email address for notifications
#SBATCH --mail-type=ALL

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
# ulimit -s 204800
ulimit -s unlimited
ulimit -c 0

# Execute OpenMP programs
srun ./ml_rttov
