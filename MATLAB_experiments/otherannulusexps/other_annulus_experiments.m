% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..


ir = 0.1;
[W,F,bi,bo] = annulus(8,ir);
V = [W zeros(size(W,1),1)];
nrefs = 6;

edgelens = nan(1,nrefs);
err_crof = nan(1,nrefs);
err_dec = nan(1,nrefs);
err_bergou = nan(1,nrefs);

for refi=1:nrefs
    
    [~,Qdec] = biharmonic_coordinates(V, F, [], 'PaperVersion',true);
    [~,Qbergou] = biharmonic_coordinates(V, F, [], 'PaperVersion',false);
    [Qcrof, L, K, M, D] = curved_hessian(V,F);
    
    z_crof = min_quad_with_fixed(Qcrof, [], [bi;bo], ...
        [ones(size(bi,1),1);zeros(size(bo,1),1)]);
    z_dec = min_quad_with_fixed(Qdec, [], [bi;bo], ...
        [ones(size(bi,1),1);zeros(size(bo,1),1)]);
    z_bergou = min_quad_with_fixed(Qbergou, [], [bi;bo], ...
        [ones(size(bi,1),1);zeros(size(bo,1),1)]);
    Power = @(x,y) x.^y;
    z_exf = @(r) (2*Power(ir,2)*log(ir).*(-1 + Power(r,2) + 2*log(r)) - ...
        (-1 + Power(ir,2)).*(3 - 3*Power(r,2) + 2*Power(r,2).*log(r)))./ ...
        (3*Power(-1 + Power(ir,2),2) + 4*Power(ir,2).*Power(log(ir),2));
    z_exact = z_exf(normrow(V(:,1:2)));
    
    M = massmatrix(V,F);
    err_crof(refi) = (z_crof-z_exact)'*M*(z_crof-z_exact);
    err_dec(refi) = (z_dec-z_exact)'*M*(z_dec-z_exact);
    err_bergou(refi) = (z_bergou-z_exact)'*M*(z_bergou-z_exact);

    E = edges(F);
    edgelens(refi) = mean(normrow(V(E(:,2),:) - V(E(:,1),:)));
    
    if refi<nrefs
        [V,F] = loop(V,F,1);
        newb = unique(outline(F));
        bi = newb(normrow(V(newb,:)) < 0.5);
        bo = newb(normrow(V(newb,:)) > 0.5);
        V(bi,:) = ir * V(bi,:)./normrow(V(bi,:));
        V(bo,:) = V(bo,:)./normrow(V(bo,:));
    end
    
end

%Plot errors
close;
figure;
set(gcf,'WindowStyle','normal');
loglog(edgelens, err_crof, '-ob', ...
    edgelens, err_dec, '-og', ...
    edgelens, err_bergou, '-oy', ...
    edgelens,edgelens, '--k');
set(gca, 'XDir','reverse');
title('Annulus convergence plot (edgelen)');
xlabel('average edge length');
ylabel('error');
legend('l2 error of CROF Hessian', ...
    'l2 error of DEC Hessian', ...
    'l2 error of Bergou Hessian', ...
    'h');

saveas(gcf, 'otherannulus-errors-edgelen.png');
saveas(gcf, 'otherannulus-errors-edgelen.eps');


names = {'exact','crof','bergou', 'dec'};
zs = {z_exact, z_crof, z_bergou, z_dec};
%Plot annulus solutions
for i=1:numel(names)
    close;
    figure;
    hold on;
    t = {};
    Vplot = V;
    Vplot(:,3) = zs{i};
    params = {'FaceLighting','gouraud', ...
        'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
        'AmbientStrength',0.6, 'EdgeColor','none'};
    t{end+1} = tsurf(F, Vplot, 'CData', zs{i}, params{:});
    if strcmp(names{i},'exactred')
        colormap(cbrewer('Reds',30));
    else
        colormap(cbrewer('YlGn',30));
    end
    shading interp;
    view(3);
    camproj('orth');
    axis equal;
    axis off;
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    %set(gcf,'Pos',[0 0 1920 1000]);
    camlight('left');
    title(['Annulus solution with ' names{i}]);
    saveas(gcf, sprintf('otherannulus-%s.png',names{i}));
end
