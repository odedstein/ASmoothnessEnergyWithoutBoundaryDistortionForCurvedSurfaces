// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "angle_defect.h"

#include <igl/squared_edge_lengths.h>
#include <igl/internal_angles.h>
#include <igl/boundary_loop.h>


template <typename DerivedV, typename DerivedF, typename DerivedK>
IGL_INLINE void
chp::angle_defect(
                  const Eigen::MatrixBase<DerivedV>& V,
                  const Eigen::MatrixBase<DerivedF>& F,
                  Eigen::PlainObjectBase<DerivedK>& K,
                  const chp::AngleDefectType type)
{
    using Scalar = typename DerivedV::Scalar;
    using Vec3 = Eigen::Matrix<Scalar, 3, 1>;
    
    using namespace igl;
    
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    l_sq;
    squared_edge_lengths(V, F, l_sq);
    
    std::vector<std::vector<typename DerivedF::Scalar> > b;
    igl::boundary_loop(DerivedF(F), b); //This DerivedF(F) is really annoying, libigl should change the function interface to MatrixBase
    
    std::map<int, Scalar> theta_b;
    for(const auto& loop : b) {
        for(int i=0; i<loop.size(); ++i) {
            const int v=loop[i], vnext=loop[i<loop.size()-1?i+1:0],
            vprev=loop[i>0?i-1:loop.size()-1];
            if(type==ANGLEDEFECTTYPE_GHOSTTRIANGLE) {
                const Vec3 vec1=V.row(vprev)-V.row(v),
                vec2=V.row(vnext)-V.row(v);
                const Vec3 cross = vec1.cross(vec2).normalized();
                theta_b[v] = atan2(vec1.cross(vec2).dot(cross), vec1.dot(vec2));
                if(theta_b[v]<0)
                    theta_b[v] += 2.*M_PI;
            } else {
                //This is irrelevant if we don't use ghost triangles
                theta_b[v] = 0;
            }
        }
    }
    
    angle_defect_using_l_sq(F, l_sq, theta_b, K, type);
}


template <typename DerivedL_sq, typename DerivedF, typename DerivedK,
typename IntTheta_B, typename ScalarTheta_B>
IGL_INLINE void
chp::angle_defect_using_l_sq(
                             const Eigen::MatrixBase<DerivedF>& F,
                             const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                             const std::map<IntTheta_B, ScalarTheta_B>& theta_b,
                             Eigen::PlainObjectBase<DerivedK>& K,
                             const chp::AngleDefectType type)
{
    using namespace igl;
    
    DerivedL_sq theta;
    internal_angles_using_squared_edge_lengths(l_sq, theta);
    angle_defect_intrinsic(F, theta, theta_b, K, type);
}


template <typename DerivedTheta, typename DerivedF, typename DerivedK,
typename IntTheta_B, typename ScalarTheta_B>
IGL_INLINE void
chp::angle_defect_intrinsic(
                            const Eigen::MatrixBase<DerivedF>& F,
                            const Eigen::MatrixBase<DerivedTheta>& theta,
                            const std::map<IntTheta_B, ScalarTheta_B>& theta_b,
                            Eigen::PlainObjectBase<DerivedK>& K,
                            const chp::AngleDefectType type)
{
    using namespace igl;
    
    K.resize(F.maxCoeff()+1, 1);
    K.setZero();
    
    //Interior angle defect
    for(int f=0; f<F.rows(); ++f)
        for(int v=0; v<3; ++v)
            K(F(f,v)) -= theta(f,v);
    
    //Interior angle defect
    //I am a bit paranoid about numerics here, so I won't just add 2pi and
    // subtract it again in a loop over the boundary later
    for(int v=0; v<K.rows(); ++v) {
        if(theta_b.count(v)==0)
            K[v] += 2.*M_PI;
    }
    
    //Boundary angle defect
    for(const auto& t : theta_b) {
        const IntTheta_B v = t.first;
        const ScalarTheta_B theta = t.second;
        switch(type) {
            case ANGLEDEFECTTYPE_GAUSSBONNET:
                K[v] += M_PI;
                break;
            case ANGLEDEFECTTYPE_ZEROATBOUNDARY:
                K[v] = 0;
                break;
            case ANGLEDEFECTTYPE_GHOSTTRIANGLE:
                K[v] += 2.*M_PI - theta;
                break;
            default:
                assert(false && "This type is not specified.");
        }
    }
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void chp::angle_defect<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, chp::AngleDefectType);
template void chp::angle_defect_using_l_sq<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, int, double>(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, std::__1::map<int, double, std::__1::less<int>, std::__1::allocator<std::__1::pair<int const, double> > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, chp::AngleDefectType);
#endif
