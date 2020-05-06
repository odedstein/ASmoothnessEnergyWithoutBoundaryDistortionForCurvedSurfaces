// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef CHP_CROF_CURVATURE_CORRECTION_H
#define CHP_CROF_CURVATURE_CORRECTION_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <Eigen/Sparse>


namespace chp
{
    // Computes the CROF curvature correction term
    //
    // Inputs:
    //  V, F: input mesh
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //  type: the type of curvature correction
    //
    // Outputs:
    //  K: computed curvature correction term
    //  E, oE: these are computed if they are not present.
    
    enum CurvatureCorrectionType {
        CURVATURECORRECTIONTYPE_ANGLEDEFECT = 0,
        CURVATURECORRECTIONTYPE_DEFAULT = 0,
        NUM_CURVATURECORRECTIONTYPE = 1
    };
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename ScalarK>
    IGL_INLINE void
    crof_curvature_correction(
                              const Eigen::MatrixBase<DerivedV>& V,
                              const Eigen::MatrixBase<DerivedF>& F,
                              const Eigen::MatrixBase<DerivedE>& E,
                              const Eigen::MatrixBase<DerivedOE>& oE,
                              Eigen::SparseMatrix<ScalarK>& K,
                              const CurvatureCorrectionType type =
                              CURVATURECORRECTIONTYPE_DEFAULT);
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename ScalarK>
    IGL_INLINE void
    crof_curvature_correction(
                              const Eigen::MatrixBase<DerivedV>& V,
                              const Eigen::MatrixBase<DerivedF>& F,
                              Eigen::PlainObjectBase<DerivedE>& E,
                              Eigen::PlainObjectBase<DerivedOE>& oE,
                              Eigen::SparseMatrix<ScalarK>& K,
                              const CurvatureCorrectionType type =
                              CURVATURECORRECTIONTYPE_DEFAULT);
    
    
    // Version that uses intrinsic quantities as input
    //
    // Inputs:
    //  F: input mesh connectivity
    //  l_sq: squared edge lengths of each halfedge
    //  theta: the tip angles at each halfedge
    //  dA: the double area of each face
    //  kappa: the gaussian curvature at each vertex
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //  type: the type of curvature used or the integration policy
    //        (if curvature given)
    //
    // Outputs:
    //  K: computed curvature correction term
    
    template <typename DerivedF, typename DerivedL_sq, typename DerivedE,
    typename DerivedOE, typename ScalarK>
    IGL_INLINE void
    crof_curvature_correction_intrinsic(
                                        const Eigen::MatrixBase<DerivedF>& F,
                                        const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                        const Eigen::MatrixBase<DerivedE>& E,
                                        const Eigen::MatrixBase<DerivedOE>& oE,
                                        Eigen::SparseMatrix<ScalarK>& K,
                                        const CurvatureCorrectionType type =
                                        CURVATURECORRECTIONTYPE_DEFAULT);
    
    enum CurvatureCorrectionIntegrationPolicy {
        CURVATURECORRECTIONINTEGRATIONPOLICY_DELTAFCTS = 0,
        CURVATURECORRECTIONINTEGRATIONPOLICY_DEFAULT = 0,
        NUM_CURVATURECORRECTIONINTEGRATIONPOLICY = 1
    };
    
    template <typename DerivedF, typename DerivedL_sq, typename DerivedKappa,
    typename DerivedE, typename DerivedOE, typename ScalarK>
    IGL_INLINE void
    crof_curvature_correction_intrinsic(
                                        const Eigen::MatrixBase<DerivedF>& F,
                                        const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                        const Eigen::MatrixBase<DerivedKappa>& kappa,
                                        const Eigen::MatrixBase<DerivedE>& E,
                                        const Eigen::MatrixBase<DerivedOE>& oE,
                                        Eigen::SparseMatrix<ScalarK>& K,
                                        const
                                        CurvatureCorrectionIntegrationPolicy
                                        type =
                                        CURVATURECORRECTIONINTEGRATIONPOLICY_DEFAULT);
    
    template <typename DerivedF, typename DerivedL_sq, typename DerivedTheta,
    typename DerivedKappa, typename DerivedE, typename DerivedOE,
    typename ScalarK>
    IGL_INLINE void
    crof_curvature_correction_deltafcts(
                                        const Eigen::MatrixBase<DerivedF>& F,
                                        const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                        const Eigen::MatrixBase<DerivedTheta>& theta,
                                        const Eigen::MatrixBase<DerivedKappa>& kappa,
                                        const Eigen::MatrixBase<DerivedE>& E,
                                        const Eigen::MatrixBase<DerivedOE>& oE,
                                        Eigen::SparseMatrix<ScalarK>& K);
    
    
}


#ifndef IGL_STATIC_LIBRARY
#  include "crof_curvature_correction.cpp"
#endif

#endif
