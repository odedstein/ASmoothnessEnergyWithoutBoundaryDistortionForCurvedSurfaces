% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh and construct exact solution
ndivs = 1;
[V,F] = readOBJ('camelhead.obj');V(:,1)=-V(:,1);
% [V,F] = readOBJ('tower.obj');V(:,1)=-V(:,1);F=flip_ears(V,F);
[V,F] = loop(V,F,ndivs);
V = V-mean(V);
V = V / norm(max(V)-min(V));
phix = pi/2; phiy = 0; phiz = 0;
R = [1 0 0; 0 cos(phix) -sin(phix); 0 sin(phix) cos(phix)] * ...
    [sin(phiy) 0 cos(phiy); 0 1 0; cos(phiy) 0 -sin(phiy)] * ...
    [cos(phiz) -sin(phiz) 0; sin(phiz) cos(phiz) 0; 0 0 1];
V = V*R';
n = size(V,1);

% Pick exact solution - camel
tvec = [1; 0.2; 0];
tvec = tvec/norm(tvec);
f = @(x) (x*tvec < -0.1)*1.1 ...
    + (x*tvec > 0.2)*0.4;
f_func = f(V);
alpha = 5e5;
planfac = 1.5;
lapfac = 1;

% Pick exact solution - tower
% f_func = (V(:,3)>0.225 | V(:,3)<-0.16 | normrow(V(:,1:2))>0.12)*1;
% alpha = 1.15e5;
% planfac = 20;
% lapfac = 0.1;


% Solve the smoothing problems
L = cotmatrix(V,F);
M = massmatrix(V,F);
Q_lap = L' * (M\L);
z_lap = min_quad_with_fixed(0.5*(Q_lap + lapfac*alpha*M), -lapfac*alpha*M*f_func, [], []);
Q_planarhess = hessian_squared(V,F);
z_planar_hess = min_quad_with_fixed(0.5*(Q_planarhess + planfac*alpha*M), -planfac*alpha*M*f_func, [], []);
Q_curvedhess = curved_hessian(V,F);
z_curved_hess = min_quad_with_fixed(0.5*(Q_curvedhess + alpha*M), -alpha*M*f_func, [], []);

%Camel
azel = [50,20];
%Tower
% azel = [110,3];

zs = {f_func, ...
    z_lap, z_planar_hess, z_curved_hess};
names = {'input', ...
    'laplacian-zn', 'planar hess', 'curved hess'
    };

for pi=1:numel(zs)
    close;
    figure;
    hold on;
    t = {};
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.2, 'DiffuseStrength',0.4, 'AmbientStrength',0.4, ...
        'EdgeColor','none'};
    t{end+1} = tsurf(F, V, 'CData', zs{pi}, params{:});
    %Camel
    l = light('Position',[3 -3 3],'Style','infinite');
    %Tower
    % l = light('Position',[1 -1 3],'Style','infinite');
    add_shadow(t,l,'Color',[0.8 0.8 0.8],'Fade','infinite','Nudge',0);

    %Tower
    %nc = 20;
    %Camel
    nc = 18;
    CM = cbrewer('OrRd'  ,nc);
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
    %set(gcf,'Pos',[0 0 2400 1200]);
    camlight('left');
    title(['Smoothing, ' names{pi}]);
    
    [bs,~] = split_backfacing(t);
    cellfun(@(b) set(b,'DiffuseStrength',0.05,'SpecularStrength',0.0,'AmbientStrength',0.2),bs);
    set(gca,'Pos',[ 0 0 1 1]);
    
    print(gcf, ['smoothing-example-' names{pi} '.png'],'-dpng','-r300');
    
end
