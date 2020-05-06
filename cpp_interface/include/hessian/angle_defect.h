// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ANGLE_DEFECT_H
#define IGL_ANGLE_DEFECT_H

#include <igl/igl_inline.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <map>


namespace chp
{
    
    enum AngleDefectType {
        ANGLEDEFECTTYPE_GAUSSBONNET = 0,
        ANGLEDEFECTTYPE_ZEROATBOUNDARY = 1,
        ANGLEDEFECTTYPE_GHOSTTRIANGLE = 2,
        ANGLEDEFECTTYPE_DEFAULT = 2,
        NUM_ANGLEDEFECTTYPE = 3
    };
    
    // Computes the CROF curvature correction term
    //
    // Inputs:
    //  V, F: input mesh
    //  type: the type of angle defect you want to compute
    //
    // Outputs:
    //  K: computed angle defect at each vertex
    
    template <typename DerivedV, typename DerivedF, typename DerivedK>
    IGL_INLINE void
    angle_defect(
                 const Eigen::MatrixBase<DerivedV>& V,
                 const Eigen::MatrixBase<DerivedF>& F,
                 Eigen::PlainObjectBase<DerivedK>& K,
                 const AngleDefectType type = ANGLEDEFECTTYPE_DEFAULT);
    
    
    // Versions using only intrinsic values
    
    // Using squared edge lenghts
    // Inputs:
    //  l_sq: squared length of each halfedge
    //  theta_b: a map that contains each boundary vertex and the associated
    //            angle of its ghost triangle (this is always needed to know
    //            what the boundary is, even if the actual entries are not used)
    
    template <typename DerivedL_sq, typename DerivedF, typename DerivedK,
    typename IntTheta_B, typename ScalarTheta_B>
    IGL_INLINE void
    angle_defect_using_l_sq(
                            const Eigen::MatrixBase<DerivedF>& F,
                            const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                            const std::map<IntTheta_B, ScalarTheta_B>& theta_b,
                            Eigen::PlainObjectBase<DerivedK>& K,
                            const AngleDefectType type
                            = ANGLEDEFECTTYPE_DEFAULT);
    
    
    // Using squared angles
    // Inputs:
    //  theta: angle at each halfedge
    //  theta_b: a map that contains each boundary vertex and the associated
    //            angle of its ghost triangle (this is always needed to know
    //            what the boundary is, even if the actual entries are not used)
    
    template <typename DerivedTheta, typename DerivedF, typename DerivedK,
    typename IntTheta_B, typename ScalarTheta_B>
    IGL_INLINE void
    angle_defect_intrinsic(
                            const Eigen::MatrixBase<DerivedF>& F,
                            const Eigen::MatrixBase<DerivedTheta>& theta,
                            const std::map<IntTheta_B, ScalarTheta_B>& theta_b,
                            Eigen::PlainObjectBase<DerivedK>& K,
                            const AngleDefectType type
                            = ANGLEDEFECTTYPE_DEFAULT);
    
    
}


#ifndef IGL_STATIC_LIBRARY
#  include "angle_defect.cpp"
#endif

#endif
