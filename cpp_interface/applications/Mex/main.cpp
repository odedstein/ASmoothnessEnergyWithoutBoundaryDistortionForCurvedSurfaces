// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifdef LIBIGL_MATLAB_ENABLED

// Std lib
#include <iostream>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Project
#include <hessian/crof_laplacian.h>
#include <hessian/crof_curvature_correction.h>
#include <hessian/crof_differential.h>
#include <hessian/crof_massmatrix.h>

// MATLAB
#include <mex.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/MexStream.h>
#include <igl/C_STR.h>


void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    igl::matlab::MexStream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);
    
    igl::matlab::mexErrMsgTxt(nlhs>=1,
                              "At least one output argument is needed.");
    igl::matlab::mexErrMsgTxt(nlhs<=7,
                              "At most seven output arguments are allowed.");
    
    igl::matlab::mexErrMsgTxt(nrhs == 2,
                              "The number of input arguments must be 2.");
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    const int dim = mxGetN(prhs[0]);
    igl::matlab::mexErrMsgTxt(dim == 3 || dim == 2,
                              "Mesh vertex list must be #V by 2 or 3 list of vertex positions");
    igl::matlab::mexErrMsgTxt(mxGetN(prhs[1]) == 3,
                              "Faces must have three vertices");
    igl::matlab::parse_rhs_double(prhs+0,V);
    igl::matlab::parse_rhs_index(prhs+1,F);
    
    
    Eigen::SparseMatrix<double> Q, L, K, M, D;
    Eigen::MatrixXi E, oE;
    
    //Construct matrices
    chp::crof_massmatrix(V, F, E, oE, M);
    chp::crof_differential(V, F, E, oE, D);
    chp::crof_laplacian(V, F, E, oE, L);
    std::vector<Eigen::Triplet<double> > tripletListMi;
    for(int k=0; k<M.outerSize(); ++k)
        for(Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
            if(it.value() > 0)
                tripletListMi.emplace_back(it.row(), it.col(), 1./it.value());
    Eigen::SparseMatrix<double> Mi(M.rows(), M.cols());
    Mi.setFromTriplets(tripletListMi.begin(), tripletListMi.end());
    chp::crof_curvature_correction(V, F, E, oE, K);
    Q = D.transpose()*Mi*(L + K)*Mi*D;
    
    //Output matrices
    igl::matlab::prepare_lhs_double(Q, plhs+0);
    if(nlhs>=2)
        igl::matlab::prepare_lhs_double(L, plhs+1);
    if(nlhs>=3)
        igl::matlab::prepare_lhs_double(K, plhs+2);
    if(nlhs>=4)
        igl::matlab::prepare_lhs_double(M, plhs+3);
    if(nlhs>=5)
        igl::matlab::prepare_lhs_double(D, plhs+4);
    if(nlhs>=6)
        igl::matlab::prepare_lhs_index(E, plhs+5);
    if(nlhs>=7)
        igl::matlab::prepare_lhs_double(oE, plhs+6);
    
    std::cout.rdbuf(outbuf);
}


#endif

