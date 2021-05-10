// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "orient_halfedges.h"

#include <igl/all_edges.h>
#include <igl/unique_simplices.h>


template <typename DerivedF, typename DerivedE, typename DerivedOE>
IGL_INLINE void
chp::orient_halfedges(
                      const Eigen::MatrixBase<DerivedF>& F,
                      Eigen::PlainObjectBase<DerivedE>& E,
                      Eigen::PlainObjectBase<DerivedOE>& oE)
{
    using namespace igl;
    assert(F.cols()==3 && "This only works for triangle meshes.");
    
    using Int = typename DerivedF::Scalar;
    
    const Int m = F.rows();
    
    DerivedE allE, EE;
    all_edges(F, allE);
    Eigen::Matrix<Int, Eigen::Dynamic, 1> IA, IC;
    unique_simplices(allE, EE, IA, IC);
    
    E.resize(m, 3);
    oE.resize(m, 3);
    for(int f=0; f<m; ++f) {
        for(int e=0; e<3; ++e) {
            const Int ind = f + m*e;
            E(f,e) = IC(ind);
            assert((EE(E(f,e),0)==allE(ind,0) || EE(E(f,e),0)==allE(ind,1)) &&
                   "Something is wrong in the edge matrix.");
            oE(f,e) = EE(E(f,e),0)==allE(ind,0) ? 1 : -1;
        }
    }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void chp::orient_halfedges<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
