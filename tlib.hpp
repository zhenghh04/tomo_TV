//
//  tlib.hpp
//  TV
//
//  Created by Hovden Group on 5/6/19.
//  Copyright © 2019 Jonathan Schwartz. All rights reserved.
//

#ifndef tlib_hpp
#define tlib_hpp

#include <stdio.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>

void tomography(Eigen::MatrixXf& recon, Eigen::MatrixXf& tiltSeries, Eigen::VectorXf& innerProduct, Eigen::SparseMatrix<float>& A, int beta);

float rmepsilonScalar(float input);

void rmepsilonVector(Eigen::VectorXf& input);

void parallelRay(int& Nray, Eigen::VectorXf angles);

#endif /* tlib_hpp */
