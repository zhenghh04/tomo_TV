#include "mpi_ctvlib.h"
#include "mpi.h"
#include "hdf5.h"
#include <iostream>
using namespace std; 

int main(int argc, char **argv) {
  int niter = 10;
  int ng = 10; 
  float beta = 0.25; 
  float beta_red = 0.985; 
  float eps = 0.019; 
  float r_max = 0.95;
  float alpha_red = 0.95; 
  float alpha = 0.2; 
  int SNR = 100; 
  bool noise = true; 
  bool save_recon = true; 
  float dPOCS, dp;
  cout << "checking.." << endl; 
  MPI_Init(&argc, &argv);
  int rank, nproc; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char fname[255] = "../Tilt_Series/256_au_sto.h5";
  char measure[255] = "../Measurement_Matrices/256_au_sto_measurement.h5";
  mpi_ctvlib tomo_obj;
  if (rank==0) cout << "load volume" << endl; 
  tomo_obj.loadVolume(fname);
  if (rank==0) cout << "load measurement matrix" << endl; 
  tomo_obj.loadMeasurementMatrix(measure);
  if (rank==0) cout << "normalization" << endl; 
  tomo_obj.normalization();
  for (int i = 0; i<niter; i++) {
    if (rank==0) cout << "iter " << i << endl; 
    tomo_obj.copy_recon();
    tomo_obj.sART(beta, -1); 
    tomo_obj.positivity();
    beta *= beta_red; 
    tomo_obj.forwardProjection(-1);
    if (i==0) {
      dPOCS = tomo_obj.matrix_2norm()*alpha; 
      dp = dPOCS/alpha; 
    } else {
      dp = tomo_obj.matrix_2norm();
    }
    float dd = tomo_obj.vector_2norm();
    float tv = tomo_obj.tv_3D();
    float rmse = tomo_obj.rmse();
    tomo_obj.copy_recon();
    tomo_obj.tv_gd_3D(ng, dPOCS);
    float dg = tomo_obj.matrix_2norm();
    if (dg > dp*r_max and dd > eps)
      dPOCS *= alpha_red; 
  }
  MPI_Finalize();
  return 0; 
}
