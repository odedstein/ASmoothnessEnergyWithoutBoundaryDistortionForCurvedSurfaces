% Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com> and
% Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.

function [Q,H,M4,D,A,G,int] = hessian_squared(V,F,varargin)
  % HESSIAN_SQUARED Construct a matrix to compute the integrated squared
  % Hessian of a function over a triangle mesh using the mixed finite element
  % method.
  % 
  % [Q] = hessian_squared(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %   extraint #any boundary vertices that are explicitly to be treated as
  %   the interior
  % Outputs:
  %   Q  #V by #V sparse matrix so that X'*Q*X measures the integrated squared
  %     Hessian energy of a scalar function X
  %
  
  extraints = [];
  if(length(varargin)>0)
      extraints = varargin{1};
  end

  % Number of faces
  m = size(F,1);
  % Number of vertices
  n = size(V,1);
  % Number of dimension
  dim = size(V,2);
  m = size(F,1);
  n = size(V,1);
  G = grad(V,F);
  assert(size(G,1) == dim*m,'Gradient should equal dim*m');
  M = massmatrix(V,F);
  M4 = repdiag(M,dim^2);
  % Block transpose of G
  GG = sparse(m,0);
  for d = 1:dim
    GG = [GG G((d-1)*m+(1:m),:)];
  end
  D = repdiag(GG,dim);

  A = repdiag(diag(sparse(doublearea(V,F)*0.5)),dim);
  H = D'*A*G;

  b = unique(outline(F));
  int = setdiff(1:n,b);
  int = [int(:); extraints];
  int4 = (0:(dim^2-1))*n + int;
  notint4 = setdiff(1:(n*dim^2), int4);
  H(notint4,:) = 0;
  Q = H' * (M4 \ H);
end
