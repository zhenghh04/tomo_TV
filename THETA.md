# This is some instruction particularly for running the code on Theta

## loading python module and install requried packages

* module load intelpython36 cray-hdf5-parallel
* pip install scipy scikit-image Pillow pybind11 --user


## Compile the C++ code
* cd 3D/Utils, modify the parameters (number of iterations) in sim_tv_mpi.C if necessary 
* cp make.inc.craype make.inc
* make 

## run the job
* aprun ... ./sim_tv -d DATA -m MEASUREMENT_MATRIX -o OUTPUT_FOR_RECONSTRUCTED_IMAGE
* An example of submission script is qsub.sc


### Running jobs on the login node
* loading the mpi module on the login node
  $ module load /soft/datascience/intelmpi/intelmpi-login

* Specify HDF5_ROOT
  
  export HDF5_ROOT=...
  
  epport LD_LIBRARY_PATH=$HDF5_ROOT/lib:$LD_LIBRARY_PATH

* build the code
  
  cp make.inc.login make.inc
     
  make

* run 
  mpirun -np ...
  