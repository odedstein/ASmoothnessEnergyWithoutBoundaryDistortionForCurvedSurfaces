// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "crof_differential.h"

#include <vector>

#include "orient_halfedges.h"

#include <igl/doublearea.h>
#include <igl/squared_edge_lengths.h>


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename ScalarD>
IGL_INLINE void
chp::crof_differential(
                       const Eigen::MatrixBase<DerivedV>& V,
                       const Eigen::MatrixBase<DerivedF>& F,
                       const Eigen::MatrixBase<DerivedE>& E,
                       const Eigen::MatrixBase<DerivedOE>& oE,
                       Eigen::SparseMatrix<ScalarD>& D)
{
    using namespace igl;
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    l_sq;
    squared_edge_lengths(V, F, l_sq);
    crof_differential_intrinsic(F, l_sq, E, oE, D);
}


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename ScalarD>
IGL_INLINE void
chp::crof_differential(
                       const Eigen::MatrixBase<DerivedV>& V,
                       const Eigen::MatrixBase<DerivedF>& F,
                       Eigen::PlainObjectBase<DerivedE>& E,
                       Eigen::PlainObjectBase<DerivedOE>& oE,
                       Eigen::SparseMatrix<ScalarD>& D)
{
    using namespace igl;
    if(E.rows()!=F.rows() || E.cols()!=F.cols() || oE.rows()!=F.rows() ||
       oE.rows()!=F.cols())
        orient_halfedges(F, E, oE);
    
    crof_differential(V, F,
                    const_cast<const Eigen::PlainObjectBase<DerivedE>& >(E),
                    const_cast<const Eigen::PlainObjectBase<DerivedOE>& >(oE),
                    D);
}


template <typename DerivedF, typename DerivedL_sq, typename DerivedE,
typename DerivedOE, typename ScalarD>
IGL_INLINE void
chp::crof_differential_intrinsic(
                                 const Eigen::MatrixBase<DerivedF>& F,
                                 const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                 const Eigen::MatrixBase<DerivedE>& E,
                                 const Eigen::MatrixBase<DerivedOE>& oE,
                                 Eigen::SparseMatrix<ScalarD>& D)
{
    using namespace igl;
    Eigen::Matrix<typename DerivedL_sq::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    dA;
    DerivedL_sq l_sqrt = l_sq.array().sqrt().matrix();
    doublearea(l_sqrt, dA);
    crof_differential_intrinsic(F, l_sq, dA, E, oE, D);
}


template <typename DerivedF, typename DerivedL_sq, typename DeriveddA,
typename DerivedE, typename DerivedOE, typename ScalarD>
IGL_INLINE void
chp::crof_differential_intrinsic(
                                 const Eigen::MatrixBase<DerivedF>& F,
                                 const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                                 const Eigen::MatrixBase<DeriveddA>& dA,
                                 const Eigen::MatrixBase<DerivedE>& E,
                                 const Eigen::MatrixBase<DerivedOE>& oE,
                                 Eigen::SparseMatrix<ScalarD>& D)
{
    using namespace igl;
    assert(F.cols()==3 && "Faces have three vertices");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    
    const int m = F.rows();
    const int n = F.maxCoeff() + 1;
    const int nE = E.maxCoeff() + 1;
    
    std::vector<Eigen::Triplet<ScalarD> > tripletList;
    tripletList.reserve(5*3*m);
    for(int f=0; f<m; ++f) {
        for(int e=0; e<3; ++e) {
            const int i=F(f, (e+1)%3), j=F(f, (e+2)%3), k=F(f, e);
            const ScalarD o=oE(f,e),
            eij=l_sq(f,e), ejk=l_sq(f,(e+1)%3), eki=l_sq(f,(e+2)%3); //These are squared quantities.
            
            tripletList.emplace_back(E(f,e), i, -o*dA(f)/(6.*eij));
            tripletList.emplace_back(E(f,e)+nE, i, -o*(eij+ejk-eki)/(12.*eij));
            tripletList.emplace_back(E(f,e), j, o*dA(f)/(6.*eij));
            tripletList.emplace_back(E(f,e)+nE, j, -o*(eij-ejk+eki)/(12.*eij));
            tripletList.emplace_back(E(f,e)+nE, k, o/6.);
        }
    }
    D.resize(2*nE, n);
    D.setFromTriplets(tripletList.begin(), tripletList.end());
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void chp::crof_differential<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::SparseMatrix<double, 0, int>&);
#endif
