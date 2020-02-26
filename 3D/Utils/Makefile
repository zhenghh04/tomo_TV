#Makefile
#!/bin/sh
include ./make.inc
CXXFLAGS += -DDEBUG #-xMIC-AVX512
PYCONFIG = ctvlib`python3-config --extension-suffix`
MPICONFIG = mpi_ctvlib`python3-config --extension-suffix`
all: ctvlib mpi_ctvlib test

ctvlib: ctvlib.cpp  ctvlib.hpp
	$(CXX) $(CXXFLAGS) $(EIGEN) $(PYBIND) $(PYBIND11) ctvlib.cpp -o $(PYCONFIG)

mpi_ctvlib: mpi_ctvlib.cpp mpi_ctvlib.hpp
	$(MPXX) $(CXXFLAGS) $(EIGEN) $(PYBIND) $(PYBIND11) mpi_ctvlib.cpp pywrap.cpp -o $(MPICONFIG) $(HDF5_LIBS) $(HPCTW)


test: sim_tv_mpi.C mpi_ctvlib.C mpi_ctvlib.h
	$(MPXX) -DDEBUG -I$(HDF5_ROOT)/include -O3 $(EIGEN) -o test.x mpi_ctvlib.C sim_tv_mpi.C $(HDF5_LIBS) $(HPCTW) $(HDF5_LIBS)

clean:
	rm -rf *.so *.o
