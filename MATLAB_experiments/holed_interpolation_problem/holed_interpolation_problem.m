% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh and construct exact solution
ndivs = 3;
[V,F] = readOBJ('pittsburghbridge.obj');
[V,F] = loop(V,F,ndivs);
tmpmean = mean(V);
V = V-tmpmean;
tmpscale = norm(max(V)-min(V));
V = V / tmpscale;
phix = pi/2; phiy = 0; phiz = 0;
R = [1 0 0; 0 cos(phix) -sin(phix); 0 sin(phix) cos(phix)] * ...
    [sin(phiy) 0 cos(phiy); 0 1 0; cos(phiy) 0 -sin(phiy)] * ...
    [cos(phiz) -sin(phiz) 0; sin(phiz) cos(phiz) 0; 0 0 1];
V = V*R';

% Load holed mesh
[VV,FF] = readOBJ('pittsburghbridge_big_holed.obj');
[VV,FF] = loop(VV,FF,ndivs);
VV = (VV-tmpmean) / tmpscale;
VV = VV*R';

% Pick points on which to interpolate the functions
pts = [...
    0.15,0.05,0; ...
    -0.23,-0.1,0; ...
    -0.04,-0.21,0; ...
    0,0.2,0; ...
    -0.03,-0.015,-0.05; ...
    0,-0.015,0.24; ...
    -0.2,0.05,0.24; ...
    ];
[b, ~, Ps] = snap_points(pts, V);
[bb, ~, ~] = snap_points(pts, VV);
bc = -Ps(:,2) + 0.5*Ps(:,1).^2;

% Solve the three interpolation problems
L = cotmatrix(V,F);
M = massmatrix(V,F);
Q_lap = L' * (M\L);
z_lap = min_quad_with_fixed(Q_lap, [], b, bc);
Q_curvedhess = curved_hessian(V,F);
z_curvedhess = min_quad_with_fixed(Q_curvedhess, [], b, bc);
LL = cotmatrix(VV,FF);
MM = massmatrix(VV,FF);
QQ_lap = LL' * (MM\LL);
zz_lap = min_quad_with_fixed(QQ_lap, [], bb, bc);
QQ_curvedhess = curved_hessian(VV,FF);
zz_curvedhess = min_quad_with_fixed(QQ_curvedhess, [], bb, bc);


% Move the Ps to plot correctly
azel = [-80,30];
PhiTheta = [0.176 -1.049];
disp = [sin(PhiTheta(2))*cos(PhiTheta(1)), ...
    sin(PhiTheta(2))*sin(PhiTheta(1)), cos(PhiTheta(2))];
Ps = Ps + 0.05*disp;

Vs = {V,V,VV,VV};
Fs = {F,F,FF,FF};
zs = {z_lap,z_curvedhess,zz_lap,zz_curvedhess};
names = {'laplacian', 'hessian', 'laplacian-holed', 'hessian-holed'};


for i=1:numel(Vs)
    close;
    figure;
    hold on;
    t = {};
    off = 1;
    t{end+1} = tsurf(Fs{i}, Vs{i}, 'CData', zs{i}, ...
        fphong, 'EdgeColor', 'none', fsoft);
    scatter3(Ps(:,1), Ps(:,2), Ps(:,3), 'SizeData', 500, ...
        'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', 'LineWidth', 3, ...
        'CData', bc);
    l = light('Position',[3 -2 -4],'Style','local');
    add_shadow(t,l,'Color',[0.8 0.8 0.8],'Fade','infinite','Nudge',0);
    nc = 30;
    CM = flipud(cbrewer('PuRd'  ,nc));
    colormap(CM);
    add_isolines(t,'LineWidth',2);
    colorbar;
    view(azel);
    camproj('persp');
    axis equal;
    axis off;
    apply_ambient_occlusion([t{:}],'AddLights',false,'SoftLighting',true,...
        'Unoriented',true);
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    camlight;
    title(['Scattered data interpolation ' names{i}]);
    
    [~,bs] = split_backfacing(t);
    cellfun(@(b) set(b,'DiffuseStrength',0.05,'SpecularStrength',0.0,'AmbientStrength',0.2),bs);
    set(gca,'Pos',[ 0 0 1 1]);
    
    saveas(gcf, sprintf('bighole-%s-res-%d.png', names{i}, ndivs));
end
