% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..
ir = 0.2;

reftypes = {'regularref'};

for rti=1:numel(reftypes)
    
    switch reftypes{rti}
        case 'regularref'
            rads = zeros(7,1);
    end
    
    err = nan(1, numel(rads));
    err_old = nan(1, numel(rads));
    ns = nan(1, numel(rads));
    edgelens = nan(1, numel(rads));
    edgebadness = nan(1, numel(rads));
    
    if strcmp(reftypes{rti}, 'regularref')
        [V,F,bi,bo] = annulus(8,ir);
    end
    
    for i=1:numel(rads)
        
        r = normrow(V);
        Power = @(x,y) x.^y;
        z_exact = (2*Power(ir,2)*log(ir).*(-1 + Power(r,2) + 2*log(r)) - ...
            (-1 + Power(ir,2)).*(3 - 3*Power(r,2) + 2*Power(r,2).*log(r)))./ ...
            (3*Power(-1 + Power(ir,2),2) + 4*Power(ir,2).*Power(log(ir),2));
        
        n = size(V,1);
        [Q, L, K, eM, D] = curved_hessian(V,F);
        Q_old = hessian_squared(V,F);
        m = size(L,1);
        z_hess = min_quad_with_fixed(Q, [], [bi;bo], ...
            [ones(size(bi,1),1);zeros(size(bo,1),1)]);
        z_hess_old = min_quad_with_fixed(Q_old, [], [bi;bo], ...
            [ones(size(bi,1),1);zeros(size(bo,1),1)]);
        
        diff = z_exact - z_hess;
        diff_old = z_exact - z_hess_old;
        
        ns(i) = size(V,1);
        E = edges(F);
        edgelens(i) = mean(normrow(V(E(:,2),:) - V(E(:,1),:)));
        err(i) = sqrt(diff' * massmatrix(V,F) * diff);
        err_old(i) = sqrt(diff_old' * massmatrix(V,F) * diff_old);
        
        l = edge_lengths(V,F);
        circumrads = l(:,1).*l(:,2).*l(:,3) ./ ...
            sqrt((l(:,1)+l(:,2)+l(:,3)).*(l(:,2)+l(:,3)-l(:,1)).*(l(:,3)+l(:,1)-l(:,2)).*(l(:,1)+l(:,2)-l(:,3)));
        inrads = 0.5*sqrt(((l(:,2)+l(:,3)-l(:,1)).*(l(:,3)+l(:,1)-l(:,2)).*(l(:,1)+l(:,2)-l(:,3))) ./ (l(:,1)+l(:,2)+l(:,3)));
        edgebadness(i) = max(circumrads ./ inrads);
        
        if i~= numel(rads) && strcmp(reftypes{rti}, 'regularref')
            [V,F] = loop(V,F,1);
            newb = unique(outline(F));
            bi = newb(normrow(V(newb,:)) < 0.5);
            bo = newb(normrow(V(newb,:)) > 0.5);
            V(bi,:) = ir * V(bi,:)./normrow(V(bi,:));
            V(bo,:) = V(bo,:)./normrow(V(bo,:));
        end
    end
    
    exps = polyfit(log(ns),log(err),1);
    deg_of_c = exps(1);
    
    close;
    figure;
    set(gcf,'WindowStyle','normal');
    loglog(ns, err, '-ob', ...
        ns, err_old, '-or', ...
        ns,ns.^(-0.5), '--k');
    title('Annulus convergence plot (vertices)');
    xlabel('number of vertices');
    ylabel('error');
    legend('l2 error of curved Hessian', ...
        'l2 error of flat Hessian', ...
        'n^{-0.5}');
    
    saveas(gcf, ['convergence_annulus_plot-' reftypes{rti} '-verts.png']);
    saveas(gcf, ['convergence_annulus_plot-' reftypes{rti} '-verts.eps']);
    
    
    close;
    figure;
    set(gcf,'WindowStyle','normal');
    loglog(edgelens, err, '-ob', ...
        edgelens, err_old, '-or', ...
        edgelens, edgelens, '--k');
    set(gca, 'XDir','reverse');
    title('Annulus convergence plot (edgelens)');
    xlabel('average edge length');
    ylabel('error');
    legend('l2 error of curved Hessian', ...
        'l2 error of flat Hessian', ...
        'h');
    
    saveas(gcf, ['convergence_annulus_plot-' reftypes{rti} '-edgelens.png']);
    saveas(gcf, ['convergence_annulus_plot-' reftypes{rti} '-edgelens.eps']);
    
    
    z = z_hess;
    clf;
    set(gcf,'WindowStyle','normal');
    hold on;
    t = {};
    off = 1;
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
        'AmbientStrength',0.2, 'EdgeColor','none'};
    t{end+1} = tsurf(F, [V z], 'CData', z, params{:});
    nc = 15;
    CM = flipud(cbrewer('Blues'  ,nc));
    colormap(CM);
    %add_isolines(t,'LineWidth',2);
    colorbar;
    view([0, 20]);
    camproj('persp');
    axis equal;
    axis off;
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    title('annulus solution');
    
    saveas(gcf, ['convergence_annulus_pic-' reftypes{rti} '.png']);
    
end
