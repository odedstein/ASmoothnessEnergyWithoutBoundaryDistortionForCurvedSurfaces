% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh
[V,F] = readOBJ('spot.obj');
[V,F] = loop(V,F,2);
phix = pi/2;
R = [1,0,0; 0,cos(phix),-sin(phix); 0,sin(phix),cos(phix)];
V = V*R';

% Construct rhs
f = @(x) sqrt(5*cos(x(:,1)) + x(:,3).^2 + exp(-x(:,2).^2) - 3.78) ...
    + 10*normrow(x - [1 1 0]);
f_rhs = f(V);


I = speye(size(V,1),size(V,1));
L = cotmatrix(V,F);
M = massmatrix(V,F);
Mi = spdiags(1./spdiags(M,0),0,size(M,1),size(M,1));
Q_lap = L' * (Mi*L);
[merk, Lc, K, Mc, D] = curved_hessian(V,F);
Mci = spdiags(1./spdiags(Mc,0),0,size(Mc,1),size(Mc,1));
Q_missingcurvature = D' * (Mci * (Lc * (Mci  *D)));
Q_curvedhess = D' * (Mci * ((Lc+K) * (Mci * D)));

z_lap = min_quad_with_fixed(0.5*Q_lap, M*f_rhs, 1, 0);
z_missingcurvature = min_quad_with_fixed(0.5*Q_missingcurvature, M*f_rhs, 1, 0);
z_curvedhess = min_quad_with_fixed(0.5*merk, M*f_rhs, 1, 0);

clf;
params = {'EdgeColor','none','FaceColor','interp','FaceLighting', ...
    'gouraud','SpecularStrength',0.1,'DiffuseStrength',0.8, ...
    'AmbientStrength',0.6};
hold on;
t = {};
off = 2;
t{end+1} = tsurf(F, V+[1*off 0 0], 'CData', z_lap, ...
    'EdgeColor', 'none', params{:});
t{end+1} = tsurf(F, V+[2*off 0 0], 'CData', z_missingcurvature, ...
    'EdgeColor', 'none', params{:});
t{end+1} = tsurf(F, V+[3*off 0 0], 'CData', z_curvedhess, ...
    'EdgeColor', 'none', params{:});
l = light('Position',[0 -10 10],'Style','infinite');
nc = 20;
CM = cbrewer('PuBu'  ,nc);
colormap(CM);
add_isolines(t,'LineWidth',2);
t{end+1} = tsurf(F, V+[0*off 0 0], 'CData', ...
    (f_rhs-min(f_rhs))/(max(f_rhs)-min(f_rhs))*(max(z_lap)-min(z_lap)) + min(z_lap), ...
    'EdgeColor', 'none', params{:});
colorbar;
add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
azel = [180-22,17];
view(azel);
%camproj('persp');
axis equal;
axis off;
shading interp;
set(gcf,'Color','w');
set(gca,'Pos',[0 0 1 1]);
set(gcf,'Pos',[0 0 1920 1000]);
camlight;
title({'Solving a biharmonic equation with different energies.', ...
    'From right to left: the RHS (scaled), the solution using the Laplacian energy, the Hessian energy without curvature correction, and with curvature correction'});

saveas(gcf, 'with_and_without_curvature_correction.png');
