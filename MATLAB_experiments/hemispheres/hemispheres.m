% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh and construct exact solution
[V,F] = readOBJ('hemisphere.obj');
[~,c] = min(normrow(V-[0 -1 0.5]));
Z = zeros(size(V,1),1);

% Pick points on which to interpolate the functions
pts = [0,-1,0.1; ...
    1,0,0.1; -1,0,0.1; ...
    0.884661,0,1.3769; -0.884661,0,1.3769; ...
    ];
[b, ~, Ps] = snap_points(pts, V);
bc = [-1;-1;-1; 1;1];

% Solve the three interpolation problems
L = cotmatrix(V,F);
M = massmatrix(V,F);
Q_lap = L' * (M\L);
z_lap = min_quad_with_fixed(Q_lap, [], b, bc);
Q_planarhess = hessian_squared(V,F);
z_planarhess = min_quad_with_fixed(Q_planarhess, [], b, bc);
Q_curvedhess = curved_hessian(V,F);
z_curvedhess = min_quad_with_fixed(Q_curvedhess, [], b, bc);


% Move the Ps to plot correctly
azel = [-22,17];
PhiTheta = [4.4 1.3];
disp = [sin(PhiTheta(2))*cos(PhiTheta(1)), ...
    sin(PhiTheta(2))*sin(PhiTheta(1)), cos(PhiTheta(2))];
Ps = Ps + 0.15*disp;

clf;
hold on;
t = {};
off = 3.5;
t{end+1} = tsurf(F, V+[0*off 0 0], 'CData', Z, ...
    fphong, 'EdgeColor', 'none', fsoft);
scatter3(Ps(:,1)+0*off, Ps(:,2), Ps(:,3), 'SizeData', 500, ...
    'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', 'LineWidth', 3, ...
    'CData', bc);
t{end+1} = tsurf(F, V+[1*off 0 0], 'CData', z_lap, ...
    fphong, 'EdgeColor', 'none', fsoft);
t{end+1} = tsurf(F, V+[2*off 0 0], 'CData', z_planarhess, ...
    fphong, 'EdgeColor', 'none', fsoft);
t{end+1} = tsurf(F, V+[3*off 0 0], 'CData', z_curvedhess, ...
    fphong, 'EdgeColor', 'none', fsoft);
l = light('Position',[0 10 10],'Style','infinite');
add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
nc = 10;
CM = cbrewer('PiYG'  ,nc);
colormap(CM);
add_isolines(t,'LineWidth',2);
colorbar;
view(azel);
camproj('persp');
axis equal;
axis off;
set(gcf,'Color','w');
set(gca,'Pos',[0 0 1 1]);
set(gcf,'Pos',[0 0 1920 1000]);
camlight;
title('Scattered data interpolation on a hemisphere. From left to right: the interpolation points, Laplacian zero Neumann, planar Hessian, curved Hessian');

saveas(gcf, 'hemispheres.png');
