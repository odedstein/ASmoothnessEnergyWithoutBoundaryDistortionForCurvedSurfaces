% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh and construct exact solution
ndivs = 2;
[V,F] = readOBJ('helmet.obj');
[V,F] = loop(V,F,ndivs);
V = V-mean(V);
V = V / norm(max(V)-min(V));
phix = pi/2; phiy = 0; phiz = 0;
R = [1 0 0; 0 cos(phix) -sin(phix); 0 sin(phix) cos(phix)] * ...
    [sin(phiy) 0 cos(phiy); 0 1 0; cos(phiy) 0 -sin(phiy)] * ...
    [cos(phiz) -sin(phiz) 0; sin(phiz) cos(phiz) 0; 0 0 1];
V = V*R';
V(:,1) = -V(:,1);

% Pick exact solution
tvec = [1; 0.2; 0];
tvec = tvec/norm(tvec);
[pt,~] = snap_points([0.2 -0.2 0.6], V);
f = @(x) exp(-(5*normrow(x-[-0.3 0.4 -0.2])).^2/2) ...
    - 0.7*heat_geodesic(V,F,pt,0.001) ...
    + 1 - x(:,3);
raw_f = f(V);
rng(6);
f_func = raw_f + 0.3*(max(raw_f)-min(raw_f))*(rand(size(raw_f)) ...
    + 0.2*sin(50*V(:,2)).*cos(60*V(:,3))) ;

% Solve the smoothing problems
L = cotmatrix(V,F);
M = massmatrix(V,F);
Q_lap = L' * (M\L);
alpha = 5e5; %5*10^(11.2);
z_lap = min_quad_with_fixed(0.5*(Q_lap + alpha*M), -alpha*M*f_func, [], []);
Q_planarhess = hessian_squared(V,F);
planfac = 1.5;
z_planar_hess = min_quad_with_fixed(0.5*(Q_planarhess + planfac*alpha*M), -planfac*alpha*M*f_func, [], []);
Q_curvedhess = curved_hessian(V,F);
z_curved_hess = min_quad_with_fixed(0.5*(Q_curvedhess + alpha*M), -alpha*M*f_func, [], []);

azel = [-130,10];

zs = {f_func, z_lap, z_planar_hess, z_curved_hess};
names = {'noisy', 'laplacian-sol', 'planarhessian-sol', ...
    'curvedhessian-sol'};

for i=1:numel(zs)
    pause(5);
    close;
    pause(5);
    figure
    hold on;
    t = {};
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.2, 'DiffuseStrength',0.4, 'AmbientStrength',0.4, ...
        'EdgeColor','none'};
    t{end+1} = tsurf(F, V, 'CData', zs{i}, params{:});
    l = light('Position',[-3 4 4],'Style','local');
    add_shadow(t,l,'Color',[0.8 0.8 0.8],'Fade','local','Nudge',0);
    nc = 30;
    CM = flipud(cbrewer('PuBuGn'  ,nc));
    colormap(CM);
    caxis([min(f_func), max(f_func)]);
    if ~strcmp(names{i}, 'noisy')
        add_isolines(t,'LineWidth',2);
    end
    colorbar;
    view(azel);
    camproj('persp');
    apply_ambient_occlusion([t{:}],'AddLights',false,'SoftLighting',true);
    axis equal;
    axis off;
    set(gcf,'Color','w');
    %set(gcf,'Pos',[0 0 2400 1200]);
    camlight;
    title({'Smoothing a noisy function', names{i}});
    [ts,bs] = split_backfacing(t);
    cellfun(@(b) set(b,'DiffuseStrength',0.05,'SpecularStrength',0.0,'AmbientStrength',0.2),ts);
    set(gca,'Pos',[ 0 0 1 1])
    
    saveas(gcf, ['dense-smoothing-' names{i} '.png']);
    
end
