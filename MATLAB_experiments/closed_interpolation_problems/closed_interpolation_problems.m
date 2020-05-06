% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh
nsubdivs = 3;
[V,F] = readOBJ('horse.obj');
bbox = norm(max(V)-min(V));
V = V / bbox;
[V,F] = loop(V,F,nsubdivs);
phix = pi/2; phiy = 2.4; phiz = 0;
R = [1 0 0; 0 cos(phix) -sin(phix); 0 sin(phix) cos(phix)] * ...
    [sin(phiy) 0 cos(phiy); 0 1 0; cos(phiy) 0 -sin(phiy)] * ...
    [cos(phiz) -sin(phiz) 0; sin(phiz) cos(phiz) 0; 0 0 1];
V = V*R';
shadV = V*[1 0 0; 0 1 0; 0 0 0] + [0 0 min(V(:,3))];

% Pick interpolation points
% Pick points on which to interpolate the functions
pts = [...
    0.26 0.05 -0.16; ...
    0.23 -0.3 -0.15; ...
    0.3 -0.3 0.0; ...
    -0.18 -0.3 0.002; ...
    -0.08 -0.3 -0.14; ...
    -0.45 0.25 0.15; ...
    ]*R';
[b, ~, Ps] = snap_points(pts, V);
bc = [...
    0; ...
    0.25; ...
    0.25; ...
    0.75; ...
    0.75; ...
    1; ...
    ];

% Solve the three interpolation problems
L = cotmatrix(V,F);
M = massmatrix(V,F);
Q_lap = L' * (M\L);
z_lap = min_quad_with_fixed(Q_lap, [], b, bc);
Q_planarhess = hessian_squared(V,F);
z_planarhess = min_quad_with_fixed(Q_planarhess, [], b, bc);
Q_curvedhess = curved_hessian(V,F);
z_curvedhess = min_quad_with_fixed(Q_curvedhess, [], b, bc);

%Compute grads
n = size(V,1);
m = size(F,1);
G = grad(V,F);
internalAng = internalangles(V,F);
anglesum = sparse(F,1,internalAng,size(V,1),1);
avgmat = sparse([F(:,1);F(:,2);F(:,3)], [1:m 1:m 1:m]', ...
    internalAng(:)./anglesum([F(:,1);F(:,2);F(:,3)]), n,m);
grad_lap = G*z_lap;
grad_lap = avgmat*normrow([grad_lap(1:m) grad_lap((m+1):(2*m)) grad_lap((2*m+1):(3*m))]);
grad_planarhess = G*z_planarhess;
grad_planarhess = avgmat*normrow([grad_planarhess(1:m) grad_planarhess((m+1):(2*m)) grad_planarhess((2*m+1):(3*m))]);
grad_curvedhess = G*z_curvedhess;
grad_curvedhess = avgmat*normrow([grad_curvedhess(1:m) grad_curvedhess((m+1):(2*m)) grad_curvedhess((2*m+1):(3*m))]);


% Move the Ps to plot correctly
azel = [0, 25];
PhiTheta = [-1.568 1.136];
disp = [sin(PhiTheta(2))*cos(PhiTheta(1)), ...
    sin(PhiTheta(2))*sin(PhiTheta(1)), cos(PhiTheta(2))];
Ps = Ps + 0.05*disp;

zs = {z_lap,z_planarhess,z_curvedhess};
zgrads = {grad_lap,grad_planarhess,grad_curvedhess};
names = {'laplacian','planar-hess','curved-hess'};

% Plot result
for i=1:numel(zs)
    close;
    figure;
    hold on;
    off = 1.5*max(V(:,1))-min(V(:,1));
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.2, 'DiffuseStrength',0.2, 'AmbientStrength',0.8, ...
        'EdgeColor','none'};
    t = {};
    s = {};
    shadcol = [0.85 0.85 0.85];
    t{end+1} = tsurf(F,V+[0*off 0 0], 'CData',zs{i}, params{:});
    scatter3(Ps(:,1)+0*off, Ps(:,2), Ps(:,3), 'SizeData', 500, ...
        'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', 'LineWidth', 3, ...
        'CData', bc);
    s{end+1} = tsurf(F,shadV+[0*off 0 0], 'CData',0*shadV(:,3), params{:}, ...
        'FaceColor', shadcol);
    view(azel);
    cs = 30;
    CM = cbrewer('YlGnBu', ceil(1.5*cs));
    colormap(CM(1:cs,:));
    colorbar;
    add_isolines(t,'LineWidth',1);
    l = light('Position',[10 -15 20],'Style','Infinite');
    hold off;
    apply_ambient_occlusion([t{:}],'AddLights',false,'SoftLighting',false,...
        'Unoriented',true);
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    set(gca,'Visible','off');
    camlight;
    axis equal;
    title({'Scattered data interpolation',names{i}});
    saveas(gcf, ['closed_interpolation_problems' names{i} '.png']);
    
    %Plot gradient
    close;
    figure;
    hold on;
    off = 1.5*max(V(:,1))-min(V(:,1));
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.05, 'DiffuseStrength',0.8, 'AmbientStrength',0.3, ...
        'EdgeColor','none'};
    t = {};
    t{end+1} = tsurf(F,V+[0*off 0 0], 'CData',zgrads{i}, params{:});
    hold off;
    view(azel);
    l = light('Position',[10 -15 20],'Style','Infinite');
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    set(gca,'Visible','off');
    caxis([min(grad_planarhess)-0.1*abs(min(grad_planarhess)), ...
        max(grad_planarhess)+0.1*abs(max(grad_planarhess))]);
    camlight;
    axis equal;
    colormap(hot(200));
    colorbar;
    title({'Scattered data interpolation grdients',names{i}});
    saveas(gcf, ['closed_interpolation_problems' names{i} '_grads.png']);
end


