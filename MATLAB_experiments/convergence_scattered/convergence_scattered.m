% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..
ir = 0.2;

nsubdivs = 9;

[V,F] = readOBJ('hemisphere.obj');
phix = -pi/2; phiy = 0; phiz = 0;
R = [1 0 0; 0 cos(phix) -sin(phix); 0 sin(phix) cos(phix)] * ...
    [sin(phiy) 0 cos(phiy); 0 1 0; cos(phiy) 0 -sin(phiy)] * ...
    [cos(phiz) -sin(phiz) 0; sin(phiz) cos(phiz) 0; 0 0 1];
V = V*R';

S = {speye(size(V,1), size(V,1))};

Pts = [0 0 1; 1 1 0.4; -1 1 0.3; 1 -1 0.3; -1 -1 0.4];
ptf = @(x) 0.1*x(:,2) + 0.2*(x(:,3)).^3 + sin(2*pi*x(:,1));

zs = {};
ns = zeros(1,nsubdivs);
edgelens = zeros(1,nsubdivs);

for i=1:nsubdivs
    [b, ~, Ps] = snap_points(Pts, V);
    bc = ptf(Ps);
    
    n = size(V,1);
    Q = curved_hessian(V,F);
    zs{i} = min_quad_with_fixed(Q, [], b, bc);
    
    ns(i) = size(V,1);
    E = edges(F);
    edgelens(i) = mean(normrow(V(E(:,2),:) - V(E(:,1),:)));
    
    if i~=nsubdivs
        [V,F,SS] = loop(V,F,1);
        S{i+1} = speye(size(V,1), size(V,1));
        for j=1:i
            S{j} = SS*S{j};
        end
        newb = unique(outline(F));
        V(newb,:) = [V(newb,1:2) ./ normrow(V(newb,1:2)) ...
            zeros(size(newb,1),1)];
        V = V./normrow(V);
    end
end

M = massmatrix(V,F);

ns = ns(1:(end-1));
edgelens = edgelens(1:(end-1));
err = nan(1,nsubdivs-1);
for i=1:(nsubdivs-1)
    z_ups = S{i}*zs{i};
    zdiff = z_ups - zs{end};
    
    err(i) = sqrt(zdiff'*M*zdiff);
end


PhiTheta = [-2.225 1.05];
disp = [sin(PhiTheta(2))*cos(PhiTheta(1)), ...
    sin(PhiTheta(2))*sin(PhiTheta(1)), cos(PhiTheta(2))];
Ps = Ps + 0.1*disp;


exps = polyfit(log(ns),log(err),1);
deg_of_c = exps(1);


%Plot in terms of verts

close;
figure;
set(gcf,'WindowStyle','normal');
loglog(ns, err, '-ob', ...
    ns,ns.^(-0.5), '--k');
title('Scattered convergence plot (verts)');
xlabel('number of vertices');
ylabel('error');
legend('l2 error of interpolated data', ...
    'n^{-0.5}');

saveas(gcf, 'convergence-scattered-plot-verts.png');
saveas(gcf, 'convergence-scattered-plot-verts.eps');


close;
figure;
set(gcf,'WindowStyle','normal');
loglog(edgelens, err, '-ob', ...
    edgelens, edgelens, '--k');
set(gca, 'XDir','reverse');
title('Scattered convergence plot (edgelens)');
xlabel('average edge length');
ylabel('error');
legend('l2 error of interpolated data', ...
    'h');

saveas(gcf, 'convergence-scattered-plot-edgelens.png');
saveas(gcf, 'convergence-scattered-plot-edgelens.eps');


clf;
set(gcf,'WindowStyle','normal');
hold on;
t = {};
off = 1;
params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
    'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
    'AmbientStrength',0.2, 'EdgeColor','none'};
scatter3(Ps(:,1)+0*off, Ps(:,2), Ps(:,3), 'SizeData', 500, ...
    'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k', 'LineWidth', 3, ...
    'CData', bc);
t{end+1} = tsurf(F, V, 'CData', zs{end}, params{:});
l = light('Position',[0 -10 10],'Style','infinite');
add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
nc = 20;
CM = flipud(cbrewer('Greens'  ,nc));
colormap(CM);
add_isolines(t,'LineWidth',2);
colorbar;
view(3);
camproj('persp');
axis equal;
axis off;
set(gcf,'Color','w');
set(gca,'Pos',[0 0 1 1]);
set(gcf,'Pos',[0 0 1920 1000]);
title('annulus solution');

saveas(gcf, 'convergence-scattered.png');
