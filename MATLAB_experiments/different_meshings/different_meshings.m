% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

meshings = {'regularmeshing', 'irregularmeshing', ...
    'asymmirregularmeshing'};

for mi = 1:numel(meshings)
    meshing = meshings{mi};
    fprintf('doing meshing %s\n', meshing);
    
    nsubdivs = 4;
    switch meshing
        case 'regularmeshing'
            [rawV,rawF] = readOBJ('regular.obj');
        case 'irregularmeshing'
            [rawV,rawF] = readOBJ('irregular.obj');
        case 'asymmirregularmeshing'
            [rawV,rawF] = readOBJ('irregular_asymm.obj');
    end
    V = rawV;
    F = rawF;
    for i = 1:nsubdivs
        [V,F] = loop(V,F,1);
        ovb = unique(outline(F));
        V(ovb,3) = 0;
    end
    
    % Pick points on which to interpolate the functions
    pts = [0,-1,0.1; ...
        1,0,0.1; -1,0,0.1; ...
        0.884661,0,1.3769; -0.884661,0,1.3769; ...
        ];
    [b, ~, Ps] = snap_points(pts, V);
    bc = [-1;-1;-1; 1;1];
    
    n = size(V,1);
    Q = curved_hessian(V,F);
    z = min_quad_with_fixed(Q, [], b, bc);
    
    
    figure;
    clf;
    set(gcf,'WindowStyle','normal');
    hold on;
    t = {};
    params = {'FaceColor','flat', 'FaceLighting','flat', ...
        'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
        'AmbientStrength',0.6, 'LineWidth',2};
    t{end+1} = tsurf(rawF, rawV, params{:});
    colormap([0.8 0.8 0.8]);
    colorbar;
    view(2);
    axis equal;
    axis off;
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    title([meshing ' wireframe']);
    saveas(gcf,sprintf('wireframes_%s.png',meshing));
    
    % Move the Ps to plot correctly
    azel = [-22,17];
    PhiTheta = [4.4 1.3];
    disp = [sin(PhiTheta(2))*cos(PhiTheta(1)), ...
        sin(PhiTheta(2))*sin(PhiTheta(1)), cos(PhiTheta(2))];
    Ps = Ps + 0.15*disp;
    
    
    figure
    clf;
    set(gcf,'WindowStyle','normal');hold on;
    t = {};
    t{end+1} = tsurf(F, V, 'CData', z, fphong, 'EdgeColor', 'none', fsoft);
    scatter3(Ps(:,1), Ps(:,2), Ps(:,3), 'SizeData', 500, ...
        'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 3, 'CData', bc);
    l = light('Position',[0 10 10],'Style','infinite');
    add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
    nc = 10;
    CM = cbrewer('PiYG'  ,nc);
    colormap(CM);
    if(exist('caxtemp'))
        caxis(caxtemp);
    else
        caxtemp = caxis;
    end
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
    title([meshing ' interpolation']);
    
    saveas(gcf, sprintf('different_meshing_%s.png',meshing));
    
end
