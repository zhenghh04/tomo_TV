
//
//  tlib.hpp
//  TV
//
//  Created by Hovden Group on 5/6/19.
//  Copyright © 2019 Jonathan Schwartz. All rights reserved.
//

#ifndef ctvlib_hpp
#define ctvlib_hpp

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

class ctvlib 
{

typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SpMat;

public: 

    // Member Variables.
    Mat *recon, *temp_recon, *tv_recon, *original_volume;
    SpMat A;
    int Nrow, Ncol, Nslice, Ny, Nz;
    Eigen::VectorXf innerProduct;
    Mat b, g;
    SpMat M; // Diagonal Weight Matrix for SIRT;
    
	// Initializes Measurement Matrix. 
	ctvlib(int Nslice, int Nray, int Nproj);

	// Initialize Experimental Projections. 
	void setTiltSeries(Mat in);
    void setOriginalVolume(Mat in, int slice);
    void create_projections();
    void poissonNoise(int SNR);

	// Constructs Measurement Matrix.
    void loadA(Eigen::Ref<Mat> pyA);
	void normalization();

	// 2D Reconstructions
	void ART(float beta, int dyn_ind);
    void SIRT(float beta, int dyn_ind);
    void positivity();
    
    // Stochastic Reconstruction
    void sART(float beta, int dyn_ind);
    std::vector<int> rand_perm(int n);
    
	//Forward Project Reconstruction for Data Tolerance Parameter. 
	void forwardProjection(int dyn_ind);

    // Acquire local copy of reconstruction.
    void copy_recon();
    
    // Measure 2-norm of projections and reconstruction.
    float matrix_2norm();
    float vector_2norm();
    float dyn_vector_2norm(int dyn_ind);
    float rmse();
    
    // Total variation
    float tv_3D();
    float original_tv_3D();
    void tv_gd_3D(int ng, float dPOCS);
    
    // Set Slices to Zero.
    void restart_recon();
    
    // Return reconstruction to python.
    Mat getRecon(int i);
    
    // Return projections to python. 
    Mat get_projections();
};

#endif /* tlib_hpp */
