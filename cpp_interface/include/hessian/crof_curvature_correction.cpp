// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "crof_curvature_correction.h"

#include <vector>

#include "orient_halfedges.h"
#include "angle_defect.h"

#include <igl/squared_edge_lengths.h>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>
#include <igl/internal_angles.h>


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename ScalarK>
IGL_INLINE void
chp::crof_curvature_correction(
                               const Eigen::MatrixBase<DerivedV>& V,
                               const Eigen::MatrixBase<DerivedF>& F,
                               const Eigen::MatrixBase<DerivedE>& E,
                               const Eigen::MatrixBase<DerivedOE>& oE,
                               Eigen::SparseMatrix<ScalarK>& K,
                               const CurvatureCorrectionType type)
{
    using namespace igl;
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    l_sq;
    squared_edge_lengths(V, F, l_sq);
    crof_curvature_correction_intrinsic(F, l_sq, E, oE, K);
}


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename ScalarK>
IGL_INLINE void
chp::crof_curvature_correction(
                               const Eigen::MatrixBase<DerivedV>& V,
                               const Eigen::MatrixBase<DerivedF>& F,
                               Eigen::PlainObjectBase<DerivedE>& E,
                               Eigen::PlainObjectBase<DerivedOE>& oE,
                               Eigen::SparseMatrix<ScalarK>& K,
                               const CurvatureCorrectionType type)
{
    using namespace igl;
    if(E.rows()!=F.rows() || E.cols()!=F.cols() || oE.rows()!=F.rows() ||
       oE.rows()!=F.cols())
        orient_halfedges(F, E, oE);
    
    crof_curvature_correction(V, F,
                              const_cast
                              <const Eigen::PlainObjectBase<DerivedE>& >(E),
                              const_cast
                              <const Eigen::PlainObjectBase<DerivedOE>& >(oE),
                              K, type);
}


template <typename DerivedF, typename DerivedL_sq, typename DerivedE,
typename DerivedOE, typename ScalarK>
IGL_INLINE void
chp::crof_curvature_correction_intrinsic(
                                         const Eigen::MatrixBase<DerivedF>& F,
                                         const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                         const Eigen::MatrixBase<DerivedE>& E,
                                         const Eigen::MatrixBase<DerivedOE>& oE,
                                         Eigen::SparseMatrix<ScalarK>& K,
                                         const CurvatureCorrectionType type)
{
    using namespace igl;
    switch(type) {
        case CURVATURECORRECTIONTYPE_ANGLEDEFECT: {
            std::vector<std::vector<typename DerivedF::Scalar> > b;
            boundary_loop(DerivedF(F), b); //This DerivedF(F) is really annoying, libigl should change the function interface to MatrixBase
            std::map<typename DerivedF::Scalar, typename DerivedL_sq::Scalar> bmap;
            for(const auto& loop : b)
                for(auto v : loop)
                    bmap[v] = 0;
            
            Eigen::Matrix<typename DerivedL_sq::Scalar,Eigen::Dynamic,1> kappa;
            angle_defect_using_l_sq(F, l_sq, bmap, kappa,
                                    ANGLEDEFECTTYPE_ZEROATBOUNDARY);
            
            crof_curvature_correction_intrinsic(F, l_sq, kappa, E, oE, K,
                                                CURVATURECORRECTIONINTEGRATIONPOLICY_DELTAFCTS);
            break;
        } default:
            assert(false && "This option is not available");
    }
}


template <typename DerivedF, typename DerivedL_sq, typename DerivedKappa,
typename DerivedE, typename DerivedOE,
typename ScalarK>
IGL_INLINE void
chp::crof_curvature_correction_intrinsic(
                                         const Eigen::MatrixBase<DerivedF>& F,
                                         const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                         const Eigen::MatrixBase<DerivedKappa>& kappa,
                                         const Eigen::MatrixBase<DerivedE>& E,
                                         const Eigen::MatrixBase<DerivedOE>& oE,
                                         Eigen::SparseMatrix<ScalarK>& K,
                                         const
                                         CurvatureCorrectionIntegrationPolicy
                                         type)
{
    using namespace igl;
    assert(F.cols()==3 && "Faces have three vertices");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    
    const int m = F.rows();
    const int nE = E.maxCoeff() + 1;
    
    switch(type) {
        case CURVATURECORRECTIONINTEGRATIONPOLICY_DELTAFCTS: {
            DerivedL_sq theta;
            internal_angles_using_squared_edge_lengths(l_sq, theta);
            crof_curvature_correction_deltafcts(F, l_sq, theta, kappa, E, oE,
                                                K);
            break;
        } default: {
            assert(false && "This option is not available");
        }
    }
}


template <typename DerivedF, typename DerivedL_sq, typename DerivedTheta,
typename DerivedKappa, typename DerivedE, typename DerivedOE,
typename ScalarK>
IGL_INLINE void
chp::crof_curvature_correction_deltafcts(
                                         const Eigen::MatrixBase<DerivedF>& F,
                                         const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                         const Eigen::MatrixBase<DerivedTheta>& theta,
                                         const Eigen::MatrixBase<DerivedKappa>& kappa,
                                         const Eigen::MatrixBase<DerivedE>& E,
                                         const Eigen::MatrixBase<DerivedOE>& oE,
                                         Eigen::SparseMatrix<ScalarK>& K)
{
    using namespace igl;
    assert(F.cols()==3 && "Faces have three vertices");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           theta.rows()==F.rows() && theta.cols()==F.cols() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    assert(kappa.rows()==F.maxCoeff()+1 &&
           "Wrong dimension in theta or kappa");
    
    const int m = F.rows();
    const int nE = E.maxCoeff() + 1;
    
    //Divide Kappa by the actual angle sum to weigh consistently.
    // Pretty sure this should be the actual angle sum here and not just 2pi.
    DerivedTheta angleSum = DerivedTheta::Zero(kappa.rows(), 1);
    for(int f=0; f<F.rows(); ++f)
        for(int j=0; j<3; ++j)
            angleSum(F(f,j)) += theta(f,j);
    const DerivedKappa scaledKappa = kappa.array() / angleSum.array();
    
    std::vector<Eigen::Triplet<ScalarK> > tripletList;
    tripletList.reserve(10*3*m);
    for(int f=0; f<m; ++f) {
        for(int e=0; e<3; ++e) {
            const ScalarK eij=l_sq(f,e), ejk=l_sq(f,(e+1)%3),
            eki=l_sq(f,(e+2)%3); //These are squared quantities.
            const ScalarK o = oE(f,e)*oE(f,(e+2)%3);
            const int i=F(f,(e+1)%3), j=F(f,(e+2)%3), k=F(f,e);
            const ScalarK ki=scaledKappa(i)*theta(f,(e+1)%3),
            kj=scaledKappa(j)*theta(f,(e+2)%3), kk=scaledKappa(k)*theta(f,e);
            
            const ScalarK costhetaidiv = (eij-ejk+eki)/(2.*eij*eki);
            const ScalarK sinthetaidiv = sqrt( (1.-pow(eij-ejk+eki,2)/
                                                (4.*eij*eki)) / (eij*eki) );
            //I'm not quite sure if it is numerically better to use the built-in
            // sin and cos functions.
            
            const ScalarK Corrijij = (ki+kj+kk) / eij;
            tripletList.emplace_back(E(f,e), E(f,e), Corrijij);
            tripletList.emplace_back(E(f,e)+nE, E(f,e)+nE, Corrijij);
            
            const ScalarK Corrijki = -o*(ki-kj-kk)*costhetaidiv;
            tripletList.emplace_back(E(f,e), E(f,(e+2)%3), Corrijki);
            tripletList.emplace_back(E(f,(e+2)%3), E(f,e), Corrijki);
            tripletList.emplace_back(E(f,e)+nE, E(f,(e+2)%3)+nE, Corrijki);
            tripletList.emplace_back(E(f,(e+2)%3)+nE, E(f,e)+nE, Corrijki);
            
            const ScalarK Corrijkiperp = o*(ki-kj-kk)*sinthetaidiv;
            tripletList.emplace_back(E(f,e), E(f,(e+2)%3)+nE, Corrijkiperp);
            tripletList.emplace_back(E(f,(e+2)%3)+nE, E(f,e), Corrijkiperp);
            tripletList.emplace_back(E(f,e)+nE, E(f,(e+2)%3), -Corrijkiperp);
            tripletList.emplace_back(E(f,(e+2)%3), E(f,e)+nE, -Corrijkiperp);
        }
    }
    
    K.resize(2*nE, 2*nE);
    K.setFromTriplets(tripletList.begin(), tripletList.end());
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void chp::crof_curvature_correction<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::SparseMatrix<double, 0, int>&, chp::CurvatureCorrectionType);
#endif

