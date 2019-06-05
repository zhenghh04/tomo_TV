//
//  tlib.cpp
//  TV
//
//  Created by Hovden Group on 5/6/19.
//  Copyright © 2019 Jonathan Schwartz. All rights reserved.
//

#include "ctvlib.hpp"
#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <iostream>

#define PI 3.14159265359

using namespace Eigen;
using namespace std;
namespace py = pybind11;

typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SpMat;

ctvlib::ctvlib(int Nslice, int Nray, int Nproj)
{
    // 
    Nx = Nray;
    Nrow = Nray*Nproj;
    Ncol = Nray*Nray;
    A.resize(Nrow,Ncol);
    innerProduct.resize(Nrow);
    b.resize(Nrow, Nx);
}

void ctvlib::setTiltSeries(Mat in)
{
    b = in;
}

void ctvlib::parallelRay(int Nray, Eigen::VectorXf angles)
{

    //Nside = Nray = y dimension of tilt series.
    int Nside = Nray;

    float pixelWidth = 1;
    float rayWidth = 1;

    //Number of projections.
    int Nproj = angles.rows();
    int idxend = 0;

    //Initialize vectors that contain matrix elements and corresponding row/column numbers.
    //Ray coordinates at 0 degrees.
    VectorXf offsets = VectorXf::LinSpaced( Nray + 1, -(Nray - 1) * 0.5, (Nray + 1) * 0.5 ) * rayWidth;
    VectorXf xgrid = VectorXf::LinSpaced( Nray + 1, -Nray * 0.5 , Nray * 0.5) * pixelWidth;
    VectorXf ygrid = VectorXf::LinSpaced( Nray + 1, - Nray * 0.5 , Nray * 0.5 ) * pixelWidth;

    //Initialize vectors that contain matrix elements and corresponding row/column numbers
    long int max_elements = 2 * Nside * Nproj * Nray;
    VectorXf rows(max_elements), cols(max_elements), vals(max_elements);

    //Loop over projection angles.
    for(int i=0; i < Nproj; i++)
    {
        float ang = angles(i) * PI / 180;
        //Points passed by rays at current angles.
        VectorXf xrayRoated = cos(ang) * offsets;
        VectorXf yrayRoated = sin(ang) * offsets;
        xrayRoated = (xrayRoated.array().abs() < 1e-8).select(0, xrayRoated);
        yrayRoated = (yrayRoated.array().abs() < 1e-8).select(0, yrayRoated);

        float a = -sin(ang);
        a = rmepsilon(a);
        float b = cos(ang);
        b = rmepsilon(b);

        //Loop rays in current projection.
        for(int j = 0; j < Nray; j++ )
        {
            //Ray: y = tx * x + intercept.
            VectorXf t_xgrid, y_xgrid, t_ygrid, x_ygrid;
            t_xgrid = (xgrid.array() - xrayRoated(j)) / a;
            y_xgrid = b * t_xgrid.array() + yrayRoated(j);
            t_ygrid = (ygrid.array() - yrayRoated(j)) / b;
            x_ygrid = a * t_ygrid.array() + xrayRoated(j);

            //Collect all points
            long tne = t_xgrid.size() + t_ygrid.size(); // tne = total number of elements
            VectorXf t_grid(tne), xx_temp(tne), yy_temp(tne);
            t_grid << t_xgrid, t_ygrid;
            xx_temp << xgrid, x_ygrid;
            yy_temp << y_xgrid, ygrid;

            //Sort the coordinates according to intersection time.
            VectorXf I = VectorXf::LinSpaced(t_grid.size(), 0, t_grid.size()-1);
            std::sort(I.data(), I.data() + I.size(), [&t_grid] (size_t i1, size_t i2) {return t_grid[i1] < t_grid[i2];} );
            VectorXf xx(I.size()), yy(I.size());
            for(int k=0; k<t_grid.size(); k++){
                xx(k) = xx_temp(I(k));
                yy(k) = yy_temp(I(k));
            }

            // Get rid of points that are outside the image grid.
            I.resize(xx.size()), I.setZero();
            float vol_boundary = Nside/2.0 * pixelWidth;
            for(int k=0; k<xx.size(); k++)
            {
                if(xx(k) >= -vol_boundary && xx(k) <= vol_boundary)
                {    if(yy(k) >= -vol_boundary && yy(k) <= vol_boundary)
                     {
                        I(k) = 1.0;
                     }

                }
            }
            removeBadElements(xx, yy, I);

            //If the ray pass through the image grid
            if (xx.size() != 0 && yy.size() != 0)
                //Get rid of double counted points.
                I.resize(xx.size()), I.setOnes();
                for(int k=0; k<xx.size()-1; k++)
                {
                    if(abs(xx(k+1) - xx(k)) <= 1e-4)
                    {    if(abs(yy(k+1) - yy(k)) <= 1e-4)
                         {
                             I(k) = 0;
                         }
                    }
                }
                removeBadElements(xx, yy, I);

                //Calculate the length within the cell.
                tne = xx.size() - 1;
                VectorXf length(tne);
                for(int k=0; k<tne; k++)
                {
                    length(k) = sqrt( pow((xx(k+1) - xx(k)),2) + pow((yy(k+1) - yy(k)),2) );
                }
                int numvals = length.size();

                //Remove the rays that are on the boundary of the box in the
                //top or to the right of the image grid
                bool check1, check2, check;
                check1 = check2 = false;
                check = true;
                if (b == 0 && abs(yrayRoated(j) - Nside/2 * pixelWidth) < 1e-15) {
                    check1 = true;
                }
                if (a == 0 && abs(xrayRoated(j) - Nside/2 * pixelWidth) < 1e-15) {
                    check2 = true;
                }
                if (check1 || check2) {
                    check = false;
                }

                //Calculate corresponding indices in measurement matrix
                if(numvals > 0 && check == true)
                {

                    tne = xx.size() - 1;
                    VectorXf midpointsX(tne);
                    VectorXf midpointsY(tne);
                    for(int k=0; k<tne; k++)
                    {
                        midpointsX(k) = 0.5 * (xx(k) + xx(k+1));
                        midpointsY(k) = 0.5 * (yy(k) + yy(k+1));
                    }

                    midpointsX = (midpointsX.array().abs() < 1e-10).select(0, midpointsX);
                    midpointsY = (midpointsY.array().abs() < 1e-10).select(0, midpointsY);

                    //Calculate the pixel index for mid points
                    VectorXf pixelIndex(tne);
                    pixelIndex = floor((Nside/2 - midpointsY.array()/pixelWidth))*Nside + floor((midpointsX.array()/pixelWidth + Nside/2));
                    pixelIndex = (pixelIndex.array() < 0 ).select(0,pixelIndex);

                    //Create the indices to store the values to the measurement matrix
                    int idxstart = idxend;
                    idxend = idxstart + numvals;
                    int idx = 0;

                    //Store row numbers, column numbers and values
                    for(int k=idxstart; k<idxend; k++)
                    {
                        rows(k) = i * Nray + j;
                        cols(k) = pixelIndex(idx);
                        vals(k) = length(idx);
                        idx = idx + 1;
                    }
                }
        }

    }

    //Truncate excess zeros.
    for(int i=0; i <idxend; i++)
    {
        A.insert(rows(i), cols(i)) = vals(i);
    }
    A.makeCompressed();
}

Mat ctvlib::recon(Eigen::Ref<Eigen::VectorXf> recon, double beta, int s, int dyn_ind)
{
    //No dynamic reconstruction, assume fully sampled batch.
    if (dyn_ind == -1)
        dyn_ind = Nrow;
    else //Calculate how many projections were sampled.
        dyn_ind *= Nx;
    
    float a;
    for(int j=0; j < dyn_ind; j++)
    {
        a = (b(s,j) - A.row(j).dot(recon)) / innerProduct(j);
        recon += A.row(j).transpose() * a * beta;
    }
    Mat foo = recon;
    foo.resize(Nx, Nx);
    return foo;
}

float ctvlib::rmepsilon(float input)
{
    if (abs(input) < 1e-10)
        input = 0;
    return input;
}

void ctvlib::removeBadElements(Eigen::VectorXf& xx, Eigen::VectorXf& yy, Eigen::VectorXf I)
{
    //Remove elements at indices where I == 0.
    VectorXf xx_temp(xx.size()), yy_temp(xx.size());
    int ind = 0;
    for(int k=0; k<xx.size();k++)
    {
        if(I(k) != 0)
        {
            xx_temp(ind) = xx(k);
            yy_temp(ind) = yy(k);
            ind++;
        }
    }
    xx.resize(ind), yy.resize(ind);
    xx = xx_temp.head(ind);
    yy = yy_temp.head(ind);
}

void ctvlib::normalization()
{
    // Row Inner Product: RIP
    for (int i = 0; i < Nrow; i++)
    {
        innerProduct(i) = A.row(i).dot(A.row(i));
    }
}

Eigen::VectorXf ctvlib::forwardProjection(Eigen::Ref<Eigen::VectorXf> recon, int dyn_ind)
{
    //No dynamic reconstruction, assume fully sampled batch.
    if (dyn_ind == -1)
        dyn_ind = Nrow;
    else //Calculate how many projections were sampled.
        dyn_ind *= Nx;
    
    VectorXf g(dyn_ind);
    for (int i = 0; i < dyn_ind; i++)
    {
        g(i) = A.row(i).dot(recon);
    }
    return g;
}


PYBIND11_MODULE(ctvlib, m)
{
    m.doc() = "C++ Scripts for TV-Tomography Reconstructions";
    py::class_<ctvlib> ctvlib(m, "ctvlib");
    ctvlib.def(py::init<int,int, int>());
    ctvlib.def("setTiltSeries", &ctvlib::setTiltSeries, "Pass the Projections to C++ Object");
    ctvlib.def("parallelRay", &ctvlib::parallelRay, "Construct Measurement Matrix");
    ctvlib.def("recon", &ctvlib::recon, "ART Tomography");
    ctvlib.def("rowInnerProduct", &ctvlib::normalization, "Calculate the Row Inner Product for Measurement Matrix");
    ctvlib.def("forwardProjection", &ctvlib::forwardProjection, "Forward Project the Reconstructions");
}

