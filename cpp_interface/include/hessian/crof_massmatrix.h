// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef CHP_CROF_MASSMATRIX_H
#define CHP_CROF_MASSMATRIX_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <Eigen/Sparse>


namespace chp
{
    // Computes the CROF mass matrix
    //
    // Inputs:
    //  V, F: input mesh
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //
    // Outputs:
    //  M: computed mass matrix
    //  E, oE: these are computed if they are not present.
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename ScalarM>
    IGL_INLINE void
    crof_massmatrix(
                    const Eigen::MatrixBase<DerivedV>& V,
                    const Eigen::MatrixBase<DerivedF>& F,
                    const Eigen::MatrixBase<DerivedE>& E,
                    const Eigen::MatrixBase<DerivedOE>& oE,
                    Eigen::SparseMatrix<ScalarM>& M);
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename ScalarM>
    IGL_INLINE void
    crof_massmatrix(
                    const Eigen::MatrixBase<DerivedV>& V,
                    const Eigen::MatrixBase<DerivedF>& F,
                    Eigen::PlainObjectBase<DerivedE>& E,
                    Eigen::PlainObjectBase<DerivedOE>& oE,
                    Eigen::SparseMatrix<ScalarM>& M);
    
    
    // Version that uses intrinsic quantities as input
    //
    // Inputs:
    //  F: input mesh connectivity
    //  l_sq: squared edge lengths of each halfedge
    //  dA: double area of each face
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //
    // Outputs:
    //  M: computed mass matrix
    
    template <typename DerivedF, typename DerivedL_sq, typename DerivedE,
    typename DerivedOE, typename ScalarM>
    IGL_INLINE void
    crof_massmatrix_intrinsic(
                              const Eigen::MatrixBase<DerivedF>& F,
                              const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                              const Eigen::MatrixBase<DerivedE>& E,
                              const Eigen::MatrixBase<DerivedOE>& oE,
                              Eigen::SparseMatrix<ScalarM>& M);
    
    template <typename DerivedF, typename DerivedL_sq, typename DeriveddA,
    typename DerivedE, typename DerivedOE, typename ScalarM>
    IGL_INLINE void
    crof_massmatrix_intrinsic(
                              const Eigen::MatrixBase<DerivedF>& F,
                              const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                              const Eigen::MatrixBase<DeriveddA>& dA,
                              const Eigen::MatrixBase<DerivedE>& E,
                              const Eigen::MatrixBase<DerivedOE>& oE,
                              Eigen::SparseMatrix<ScalarM>& M);
    
    
}


#ifndef IGL_STATIC_LIBRARY
#  include "crof_massmatrix.cpp"
#endif

#endif
