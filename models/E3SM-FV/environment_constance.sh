source /etc/profile.d/modules.sh
module purge
module load perl/5.20.0 cmake/3.3.0 python/2.7.8 intel/15.0.1 mkl/15.0.1 mvapich2/2.1 mvapich2/2.1 pnetcdf/1.6.1 netcdf/4.3.2
export OMP_STACKSIZE=64M
export NETCDF_HOME=/share/apps/netcdf/4.3.2/intel/15.0.1/lib/../
export MKL_PATH=/share/apps/intel/2015u1/mkl//lib/intel64
export PNETCDF_PATH=/share/apps/pnetcdf/1.6.1/intel/15.0.1/mvapich2/2.1
