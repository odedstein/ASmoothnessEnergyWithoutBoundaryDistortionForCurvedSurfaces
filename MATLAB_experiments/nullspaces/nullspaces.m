% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

nspmeshes = {'cheeseman'};
nvals = 6;

for i=1:numel(nspmeshes)
    [V,F] = readOBJ([nspmeshes{i} '.obj']);
    V = (V-mean(V)) / max(max(V)-min(V));
    F = flip_ears(V,F);
    
    Qcrof = curved_hessian(V,F);
    M = massmatrix(V,F);
    [vecs_crof, vals_crof] = eigs(Qcrof, M, nvals, 'sm');
    vecs_crof = real(vecs_crof);
    vals_crof = real(vals_crof);
    
    %Eval plot
    close;
    figure;
    plot(1:nvals, diag(vals_crof), '-ob');
    title('The smallest eigenvalues of the CROF Hessian matrix');
    xlabel('which evalue');
    saveas(gcf,['nullspaces-' nspmeshes{i} '-crofvals.png']);
    saveas(gcf,['nullspaces-' nspmeshes{i} '-crofvals.eps']);
    
    %Plot eigenvectors
    for veci=1:nvals
        close;
        figure;
        hold on;
        t = {};
        params = {'FaceLighting','gouraud', ...
            'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
            'AmbientStrength',0.6, 'EdgeColor','none'};
        t{end+1} = tsurf(F, V, 'CData', vecs_crof(:,veci), params{:});
        colormap(cbrewer('Purples',15));
        shading interp;
        add_isolines(t,'LineWidth',2);
        view(2);
        camproj('orth');
        axis equal;
        axis off;
        set(gcf,'Color','w');
        set(gca,'Pos',[0 0 1 1]);
        %set(gcf,'Pos',[0 0 1920 1000]);
        title(sprintf('%s %dth eigenfunction. From left to right: DEC, CROF, planar.',nspmeshes{i},veci));
        saveas(gcf, sprintf('nullspaces-%s-eigenfunc-%d.png',nspmeshes{i},veci));
    end
end