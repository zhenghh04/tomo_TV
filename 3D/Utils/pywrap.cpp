#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "mpi_ctvlib.hpp"
namespace py = pybind11;

PYBIND11_MODULE(mpi_ctvlib, m)
{
  m.doc() = "C++ Scripts for TV-Tomography Reconstructions with OpenMPI Support";
  py::class_<mpi_ctvlib> mpi_ctvlib(m, "mpi_ctvlib");
  mpi_ctvlib.def(py::init<int,int, int>());
  mpi_ctvlib.def("get_Nslice_loc", &mpi_ctvlib::get_Nslice_loc, "Get the size of local volume");
  mpi_ctvlib.def("get_first_slice", &mpi_ctvlib::get_first_slice, "Get first slice location");
  mpi_ctvlib.def("get_rank", &mpi_ctvlib::get_rank, "Get rank id");
  mpi_ctvlib.def("get_nproc", &mpi_ctvlib::get_nproc, "Get number of processor in current communicator");
  mpi_ctvlib.def("setTiltSeries", &mpi_ctvlib::setTiltSeries, "Pass the Projections to C++ Object");
  mpi_ctvlib.def("setOriginalVolume", &mpi_ctvlib::setOriginalVolume, "Pass the Volume to C++ Object");
  mpi_ctvlib.def("create_projections", &mpi_ctvlib::create_projections, "Create Projections from Volume");
  mpi_ctvlib.def("getRecon", &mpi_ctvlib::getRecon, "Return the Reconstruction to Python");
  mpi_ctvlib.def("gather_recon", &mpi_ctvlib::gather_recon, "gather reconstruction matrix");
  mpi_ctvlib.def("save_recon", &mpi_ctvlib::save_recon, "save reconstruction matrix");
  mpi_ctvlib.def("mpi_finalize", &mpi_ctvlib::mpi_finalize, "Finalize the communicator");
  mpi_ctvlib.def("ART", &mpi_ctvlib::ART, "ART Reconstruction");
  mpi_ctvlib.def("sART", &mpi_ctvlib::sART, "Stochastic ART Reconstruction");
  mpi_ctvlib.def("SIRT", &mpi_ctvlib::SIRT, "SIRT Reconstruction");
  mpi_ctvlib.def("rowInnerProduct", &mpi_ctvlib::normalization, "Calculate the Row Inner Product for Measurement Matrix");
  mpi_ctvlib.def("positivity", &mpi_ctvlib::positivity, "Remove Negative Elements");
  mpi_ctvlib.def("set_background", &mpi_ctvlib::set_background, "Set background to be certain value");
  mpi_ctvlib.def("forwardProjection", &mpi_ctvlib::forwardProjection, "Forward Projection");
  mpi_ctvlib.def("load_A", &mpi_ctvlib::loadA, "Load Measurement Matrix Created By Python");
  mpi_ctvlib.def("copy_recon", &mpi_ctvlib::copy_recon, "Copy the reconstruction");
  mpi_ctvlib.def("matrix_2norm", &mpi_ctvlib::matrix_2norm, "Calculate L2-Norm of Reconstruction");
  mpi_ctvlib.def("vector_2norm", &mpi_ctvlib::vector_2norm, "Calculate L2-Norm of Projection (aka Vectors)");
  mpi_ctvlib.def("dyn_vector_2norm", &mpi_ctvlib::dyn_vector_2norm, "Calculate L2-Norm of Partially Sampled Projections (aka Vectors)");
  mpi_ctvlib.def("rmse", &mpi_ctvlib::rmse, "Calculate reconstruction's RMSE");
  mpi_ctvlib.def("tv", &mpi_ctvlib::tv_3D, "Measure 3D TV");
  mpi_ctvlib.def("original_tv", &mpi_ctvlib::original_tv_3D, "Measure original TV");
  mpi_ctvlib.def("tv_gd", &mpi_ctvlib::tv_gd_3D, "3D TV Gradient Descent");
  mpi_ctvlib.def("get_projections", &mpi_ctvlib::get_projections, "Return the projection matrix to python");
  mpi_ctvlib.def("poissonNoise", &mpi_ctvlib::poissonNoise, "Add Poisson Noise to Projections");
  mpi_ctvlib.def("lip", &mpi_ctvlib::lipschits, "Add Poisson Noise to Projections");
  mpi_ctvlib.def("restart_recon", &mpi_ctvlib::restart_recon, "Set all the Slices Equal to Zero");
}

