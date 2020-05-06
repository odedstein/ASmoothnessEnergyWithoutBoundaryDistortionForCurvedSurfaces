% function [Q, L, K, M, D, E, oE] = curved_hessian(V,F)
% CURVED_HESSIAN Compute the sum of solid angles of a triangle (tetrahedron)
% described by points (vectors) V
%
% [W] = winding_number(V,F,O)
% [W] = winding_number(V,F,O,'ParameterName',ParameterValue, ...)
%
% Inputs:
%  V  #V by 3 list of vertex positions
%  F  #F by 3 list of triangle indices
% Outputs:
%  Q  #V by #V sparse matrix, the CROF Hessian
%  L  2*#edges by 2*#edges sparse matrix, the CROF vector Dirichlet energy
%  K  2*#edges by 2*#edges sparse matrix, the CROF integrated scalar curvature
%      matrix
%  M  2*#edges by 2*#edges sparse matrix, the CROF mass matrix
%  D  2*#edges by #V sparse matrix, the differential matrix
%  E  #edges by 2 matrix, the edges of the mesh used for the CROF matrix
%  oE #edges by 1 matrix, the edge orientations
%

