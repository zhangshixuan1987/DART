# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
. /etc/profile.d/modules.sh
module purge 
module load cmake/3.19.6 gcc/8.1.0 intel/19.0.5 intelmpi/2019u4 netcdf/4.6.3 hdf5/1.10.5 pnetcdf/1.9.0 mkl/2019u5
export NETCDF_HOME=/share/apps/netcdf/4.6.3/intel/19.0.5/
export MKL_PATH=/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/mkl
export LD_LIBRARY_PATH=/share/apps/gcc/8.1.0/lib:/share/apps/gcc/8.1.0/lib64:/share/apps/gcc/8.1.0/lib:/share/apps/gcc/8.1.0/lib64:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64_lin/gcc4.7:/share/apps/pnetcdf/1.9.0/intel/19.0.5/intelmpi/2019u4/lib:/share/apps/netcdf/4.6.3/intel/19.0.5/lib:/share/apps/hdf5/1.10.5/serial/lib:/share/apps/intel/2019u4/compilers_and_libraries_2019.4.243/linux/mpi/intel64/libfabric/lib:/share/apps/intel/2019u4/compilers_and_libraries_2019.4.243/linux/mpi/intel64/lib/release:/share/apps/intel/2019u4/compilers_and_libraries_2019.4.243/linux/mpi/intel64/lib:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin:/share/apps/intel/2019u5/comepilers_and_libraries_2019.5.281/linux/mpi/intel64/libfabric/lib:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib/release:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/ipp/lib/intel64:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.7:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.7:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/daal/lib/intel64_lin:/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/daal/../tbb/lib/intel64_lin/gcc4.4
export I_MPI_ADJUST_ALLREDUCE=1