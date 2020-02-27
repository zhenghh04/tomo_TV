#Makefile
#!/bin/sh
include ./make.inc
CXXFLAGS += -DDEBUG
PYCONFIG = ctvlib`python3-config --extension-suffix`
MPICONFIG = mpi_ctvlib`python3-config --extension-suffix`
%.o : %.C
	$(MPXX) $(CXXFLAGS) $(EIGEN) -c $< -o $@


%.o : %.cpp
	$(MPXX) $(CXXFLAGS) $(EIGEN) -c $< -o $@


all: sim_tv ctvlib mpi_ctvlib

ctvlib: ctvlib.o  ctvlib.hpp
	$(CXX) $(CXXFLAGS) $(EIGEN) $(PYBIND) $(PYBIND11) ctvlib.o -o $(PYCONFIG)

mpi_ctvlib: mpi_ctvlib.cpp mpi_ctvlib.hpp
	$(MPXX) $(CXXFLAGS) $(EIGEN) $(PYBIND) $(PYBIND11) mpi_ctvlib.cpp -o $(MPICONFIG) $(HDF5_LIBS) $(HPCTW)

sim_tv: sim_tv_mpi.o mpi_ctvlib.o mpi_ctvlib.h
	$(MPXX) $(CXXFLAGS) $(EIGEN) -o ../sim_tv mpi_ctvlib.o sim_tv_mpi.o $(HDF5_LIBS) $(HPCTW) $(HDF5_LIBS)

clean:
	rm -rf *.so *.o *.x sim_tv_mpi ../sim_tv_mpi
