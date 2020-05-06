% Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
%
% This Source Code Form is subject to the terms of the Mozilla Public License
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.
addpath ../../cpp_interface/build/applications/Mex
addpath ..

nsubdivs = 3;
[V,F] = load_mesh('nefertiti.obj');
R = axisangle2matrix([ 0 0 1],-pi*0.1);
V = V*R;
[V,F] = loop(V,F,nsubdivs);

M = massmatrix(V,F);
L = cotmatrix(V,F);
Q = hessian_squared(V,F);
Qn = curved_hessian(V,F);

[EVh,EDh] = eigs(Q,M,4,'sm');
[EVl,EDl] = eigs(L*(M\L),M,4,'sm');
[EVc,EDc] = eigs(Qn,M,4,'sm');

t = {};
clf;
hold on;
xoff = 300;
EVl(:,2) = EVl(:,2)*sign(EVl(:,2)'*EVh(:,2));
EVc(:,2) = EVc(:,2)*sign(EVc(:,2)'*EVh(:,2));
t{end+1} = tsurf(F,V+0*[xoff 0 0],'CData',EVl(:,2),'EdgeColor','none',fphong,fsoft);
t{end+1} = tsurf(F,V+1*[xoff 0 0],'CData',EVh(:,2),'EdgeColor','none',fphong,fsoft);
t{end+1} = tsurf(F,V+2*[xoff 0 0],'CData',EVc(:,2),'EdgeColor','none',fphong,fsoft);
hold off;
axis equal;
camproj('persp');
view(0,10);

% Fred Wilson, "Grey Area (Brown Version)"
CM = [213 190 166;199 171 150; 149 126 112; 118 99 87; 72 66 62]/255;
CM = interp1(linspace(0,1,size(CM,1)),CM,linspace(0,1,20),'pchip');
colormap(CM);
caxis([-0.0029017 0.0017535]);
caxis manual;
add_isolines(t,'LineWidth',2);
l = light('Position',[2 -13 10],'Style','Infinite');
add_shadow(t,l,'Color',0.8*[1 1 1],'Fade','infinite');
camlight;
apply_ambient_occlusion([t{:}],'AddLights',false,'SoftLighting',false,...
    'Unoriented',true);
set(gca,'Visible','off');
set(gca,'pos',[0 0 1 1]);
set(gcf,'Color','w');
title('Eigenvectors. From left to right: Laplacian, planar Hessian, curved Hessian.');
saveas(gcf, 'eigenvectors.png');
