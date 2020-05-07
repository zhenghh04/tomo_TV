#include "mpi_ctvlib.h"
#include "hdf5.h"
#include <iostream>
#include "timing.h"
using namespace std; 

int main(int argc, char **argv) {
  // Number of iteration
  Timing tt;
  int niter = 10;
  // iteration in the TV loop
  int ng = 10;
  // parameter in ART reconstruction
  float beta = 0.25;
  // ART reduction 
  float beta_red = 0.985;
  // Data Tolerrance Parameter
  float eps = 0.019;
  //Reduction Criteria
  float r_max = 0.95;
  float alpha_red = 0.95; 
  float alpha = 0.2;

  int SNR = 100;
  // Outcomes
  bool noise = true; 
  bool save_recon = true; 
  float dPOCS, dp;
  int i=0;
  char fname[255] = "../Tilt_Series/256_au_sto.h5";
  char measure[255] = "../256_au_sto_measurement.h5";
  char output[255] = "../output.h5";

  mpi_ctvlib tomo_obj(&argc, &argv);
  while(i<argc) {
    if (strcmp(argv[i], "-d")==0) {
      strcpy(fname, argv[i+1]);
      i = i+2;
    } else if (strcmp(argv[i], "-m")==0) {
      strcpy(measure, argv[i+1]);
      i=i+2;
    } else if (strcmp(argv[i], "-o")==0) {
      strcpy(output, argv[i+1]);
      i=i+2;
    } else if (strcmp(argv[i], "-n")==0) {
      niter = int(atof(argv[i+1]));
      i=i+2;
    } else {
      i=i+1;
    }
  }
  int verbose = tomo_obj.get_verbose();
  if (verbose==1) {
    cout << "# ===== Command Line Inputs" << endl; 
    cout << "* Data: " << fname << endl;
    cout << "* Measurement: " << measure << endl;
    cout << "* Output: " << output << endl;
    cout << "* Niter: " << niter << endl; 
  }

  if (verbose==1) cout << "# ==== load volume ==== " << endl; 
  tomo_obj.loadVolume(fname);
  if (verbose==1) cout << "# ==== load measurement matrix" << endl; 
  tomo_obj.loadMeasurementMatrix(measure);
  if (verbose==1) cout << "# ==== normalization" << endl; 
  tomo_obj.normalization();
  if (noise)
    tomo_obj.set_background(1.0);

  if (verbose==1) cout << "# ==== create projections" << endl; 
  tomo_obj.create_projections();
  if (noise)
    tomo_obj.poissonNoise(SNR);

  if (verbose==1) cout << "# ==== original_tv_3D" << endl; 
  tomo_obj.original_tv_3D();
  for (int i = 0; i<niter; i++) {
    if (verbose==1) cout << "iter " << i << endl; 
    tomo_obj.copy_recon();

    tt.start_clock("sART");
    tomo_obj.sART(beta, -1); 
    tt.stop_clock("sART");

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
    tt.start_clock("tv_3D");
    float tv = tomo_obj.tv_3D();
    tt.stop_clock("tv_3D");
    float rmse = tomo_obj.rmse();
    tomo_obj.copy_recon();
    // TV minimization
    tt.start_clock("tv_gd_3D");
    tomo_obj.tv_gd_3D(ng, dPOCS);
    tt.stop_clock("tv_gd_3D");

    float dg = tomo_obj.matrix_2norm();
    if ((dg > dp * r_max) && (dd > eps))
      dPOCS *= alpha_red;
    if (verbose==1) {
      cout << "  - dd: " << dd  << endl;
      cout << "  - rmse: " << rmse  << endl;
      cout << "  - tv: " << tv  << endl; 
    }
  }
  tomo_obj.save_recon(output, 0);
  if (verbose==1)
    tt.PrintTiming();
  return 0; 
}
