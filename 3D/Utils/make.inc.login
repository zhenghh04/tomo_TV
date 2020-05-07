CXX = g++
MPXX = mpicxx

HDF5_INC = -I$(HDF5_ROOT)/include 
HDF5_LIBS= -L$(HDF5_ROOT)/lib -lhdf5 -lhdf5_hl 
CPATH=""
CXXFLAGS_FPIC = -O3 -fPIC -shared -fopenmp -std=c++11  -I$(HOME)/.local/include/python3.6m/ `python3 -m pybind11 --includes` $(HDF5_INC) -g
CXXFLAGS = -O3 -fopenmp -std=c++11  -I$(HOME)/.local/include/python3.6m/ `python3 -m pybind11 --includes` $(HDF5_INC) -g -DEIGEN_USE_MKL_ALL -I$(MKLROOT)/include 
EIGEN = -I/soft/libraries/eigen3/3.3.5/include/eigen3
#HPCTW = -L/soft/datascience/hpctw -lhpmprof -lbfd -liberty -lz
MKL= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl