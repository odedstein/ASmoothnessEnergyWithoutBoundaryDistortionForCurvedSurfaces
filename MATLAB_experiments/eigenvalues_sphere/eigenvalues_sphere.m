% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

%CAREFUL: This version of the script artificially increases the error for
% the 0-th eigenvalue to 1e-6 to filter out numerical noise

nsubdivs = 9;

whicheig = 5;
exacteval = [0, 2, 2, 2, 6].^2;

ns = nan(1, nsubdivs);
edgelens = nan(1, nsubdivs);

errs = {};
for i=1:whicheig
    errs{i} = nan(1, nsubdivs);
end

for i=1:nsubdivs
    
    [V,F] = subdivided_sphere(i-1);
    
    n = size(V,1);
    Q = curved_hessian(V,F);
    M = massmatrix(V,F,'full');
    
    [vec,val] = eigs(Q, M, whicheig, 'sm');
    
    ns(i) = size(V,1);
    E = edges(F);
    edgelens(i) = mean(normrow(V(E(:,2),:) - V(E(:,1),:)));
    for j=1:whicheig
        err = errs{j};
        err(i) = abs(val(j,j) - exacteval(j));
        errs{j} = err;
    end
end

%PURELY FOR PLOT CONVENIENCE: round up errs{1} to 1e-6
errs{1} = max(errs{1},1e-6);

%Plot in terms of vertex
close
figure;
set(gcf,'WindowStyle','normal');
loglog(ns, errs{1}, '-ob', ...
    ns, errs{2}, '-or', ...
    ns, errs{3}, '-og', ...
    ns, errs{4}, '-oy', ...
    ns, errs{5}, '-om', ...
    ns,ns.^(-0.5), '--k');
title('Sphere eigenvalue convergence plot (vertices)');
xlabel('number of vertices');
ylabel('error');
legend('1st eigenvalue', ...
    '2nd eigenvalue', ...
    '3rd eigenvalue', ...
    '4th eigenvalue', ...
    '5th eigenvalue', ...
    'n^{-0.5}');

saveas(gcf, 'eigenvalues_sphere_plot-verts.png');
saveas(gcf, 'eigenvalues_sphere_plot-verts.eps');


%Plot in terms of average egde length
close
figure;
set(gcf,'WindowStyle','normal');
loglog(edgelens, errs{1}, '-ob', ...
    edgelens, errs{2}, '-or', ...
    edgelens, errs{3}, '-og', ...
    edgelens, errs{4}, '-oy', ...
    edgelens, errs{5}, '-om', ...
    edgelens,edgelens, '--k');
set(gca, 'XDir','reverse');
title('Sphere eigenvalue convergence plot (edgelens), l2 error');
xlabel('average edge length');
ylabel('error');
legend('1st eigenvalue', ...
    '2nd eigenvalue', ...
    '3rd eigenvalue', ...
    '4th eigenvalue', ...
    '5th eigenvalue', ...
    'h');

saveas(gcf, 'eigenvalues_sphere_plot-edgelens.png');
saveas(gcf, 'eigenvalues_sphere_plot-edgelens.eps');


z = vec(:,whicheig);
clf;
set(gcf,'WindowStyle','normal');
hold on;
t = {};
off = 1;
params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
    'SpecularStrength',0.5, 'DiffuseStrength',0.3, ...
    'AmbientStrength',0.2, 'EdgeColor','none'};
t{end+1} = tsurf(F, V, 'CData', z, params{:});
nc = 300;
CM = flipud(cbrewer('Oranges'  ,nc));
colormap(CM);
%add_isolines(t,'LineWidth',2);
colorbar;
view(3);
camproj('persp');
axis equal;
axis off;
set(gcf,'Color','w');
set(gca,'Pos',[0 0 1 1]);
set(gcf,'Pos',[0 0 1920 1000]);
title('annulus solution');

saveas(gcf, 'eigenvalues_sphere.png');
