% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

% Load mesh and construct exact solution
nsubdivs = 2;
[V,F] = readOBJ('plane_holes.obj');
V(:,1) = -V(:,1);
[V,F] = loop(V,F,nsubdivs);

% Pick points on which to interpolate the functions
pts = [1.2,7,2; ...
    1.2,-7,1; ...
    12,0,2; ...
    -12,0,2];
rotang = 0.5;
[b, ~, Ps] = snap_points(pts, V);
bc = [1 1 -1 -1]';

% Solve the three interpolation problems
L = cotmatrix(V,F);
M = massmatrix(V,F);
Q_lap = L' * (M\L);
z_lap = min_quad_with_fixed(Q_lap, [], b, bc);
Q_planarhess = hessian_squared(V,F);
z_planarhess = min_quad_with_fixed(Q_planarhess, [], b, bc);
Q_curvedhess = curved_hessian(V,F);
z_curvedhess = min_quad_with_fixed(Q_curvedhess, [], b, bc);


azel = [-30,35];
zs = {0*ones(size(V,1),1), z_lap, z_planarhess, z_curvedhess};
names = {'initial', 'lap-z-neu', 'planar-hess', 'curved-hess'};

for ni=1:numel(zs)
    close;
    figure;
    hold on;
    t = {};
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'EdgeColor','none', 'SpecularStrength',0.3, ...
        'DiffuseStrength',0.1, 'AmbientStrength',0.7};
    t{end+1} = tsurf(F, V, 'CData', sign(zs{ni}).*(abs(zs{ni})).^1, ...
        params{:});
    %if strcmp(names{ni},'initial')
        scatter3(Ps(:,1), Ps(:,2), Ps(:,3), 'SizeData', 100, ...
            'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', ...
            'LineWidth', 3, 'CData', sign(bc).*(abs(bc)).^1);
    %end
    hold off;
    l = light('Position',10*[-1 -1 3],'Style','local');
    add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
    nc = 40;
    CM = cbrewer('RdYlBu'  ,nc);
    colormap(CM);
    caxis([-1.3,1.3]);
    add_isolines(t,'LineWidth',2);
    colorbar;
    view(azel);
    camproj('persp');
    axis equal;
    axis off;
    apply_ambient_occlusion([t{:}],'AddLights',false,'SoftLighting',false,...
        'Unoriented',true);
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    %set(gcf,'Pos',[0 0 1920 1000]);
    camlight('right');
    title(['Scattered data interpolation, ' names{ni} '.']);
    
    [~,bs] = split_backfacing(t);
    cellfun(@(b) set(b,'DiffuseStrength',0.05,'SpecularStrength',0.0,'AmbientStrength',0.2),bs);
    set(gca,'Pos',[ 0 0 1 1]);
    
    print(gcf, ['teaser_' names{ni} '.png'],'-dpng','-r300');
    
end
