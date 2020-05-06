addpath ../../../cpp/build/applications/Mex
addpath ..

convexperiments = {'blob','omega','arm','spiral','churchdome','octopus'};
%convexperiments = {'omega'};

all_edgelengths = {};
minel = 1e10;
maxel = 0;
all_ns = {};
minn = 1e10;
maxn = 0;
all_errs = {};

for convexp = 1:numel(convexperiments)
    convexperiment = convexperiments{convexp};
    
    switch convexperiment
        case 'blob'
            [W,F] = bwmesh('blob.png', 'Tol',6);
            nrefs = 7;
            zfunc = @(x) 1-x(:,1).^2-x(:,2).^2;
            f = @(x) sin(0.5*x(:,1)) - cos(0.5*x(:,2));
        case 'omega'
            [W,F] = bwmesh('omega-domain.png', 'Tol',8);
            W = W-mean(W);
            W = W*[cos(1.55) -sin(1.55); sin(1.55) cos(1.55)];
            nrefs = 7;
            zfunc = @(x) 0.08*atan2(x(:,2),x(:,1));
            f = @(x) (x(:,1).^2 + x(:,2).^2).^2;
        case 'spiral'
            [W,F] = bwmesh('spiral.png', 'Tol',10);
            W = W-mean(W);
            nrefs = 7;
            arm = max(max(W)-min(W));
            zfunc = @(x) 0.05*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi + (abs(x(:,1))+abs(x(:,2))>310).*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi<1).*2*pi + (x(:,2)>210).*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi + (abs(x(:,1))+abs(x(:,2))>325).*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi<1).*2*pi<3)*2*pi + (atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi + (abs(x(:,1))+abs(x(:,2))>325).*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi<1).*2*pi + (x(:,2)>210).*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi + (abs(x(:,1))+abs(x(:,2))>310).*(atan2(x(:,2),x(:,1)) + (atan2(x(:,2),x(:,1))<-2).*2*pi<1).*2*pi<3)*2*pi < 4).*(normrow(x)>420).*(x(:,2)>0).*2*pi);
            zfuncnopl = @(x) 0.05*(atan2(x(:,2),x(:,1)));
            zfunc = @(x) zfunc(arm*x);
            f = @(x) tan(x(:,1).*x(:,2));
        case 'arm'
            [W,F] = bwmesh('arm.png', 'Tol',8);
            nrefs = 7;
            zfunc = @(x) exp(-(x(:,1).^2+x(:,2).^2)/2);
            f = @(x) x(:,1).*x(:,2);
        case 'churchdome'
            [W,F] = bwmesh('disk.png', 'Tol',7);
            W = W-mean(W);
            nrefs = 7;
            arm = max(max(W)-min(W));
            zfunc = @(x) 0.0005*(500*exp(-0.0005*normrow(x).^2) + 800*sqrt(0.001*(360-normrow(x))));
            zfunc = @(x) zfunc(arm*x);
            f = @(x) 0.001*log(1 + x(:,1) + x(:,2));
        case 'octopus'
            [W,F] = bwmesh('octopus.png', 'Tol',6);
            W = W-mean(W);
            nrefs = 7;
            arm = max(max(W)-min(W));
            zfunc = @(x) 0.00000035*normrow(x).^2;
            zfunc = @(x) zfunc(arm*x);
            f = @(x) 0.7*x(:,1).^3 - x(:,2).^3;
    end
    W = (W-mean(W)) / max(max(W)-min(W));
    
    vals = nan(1,nrefs);
    ns = nan(1,nrefs);
    edgelens = nan(1,nrefs);
    
    for i=1:nrefs
        fprintf('computing ref %d for experiment %s\n', i, ...
            convexperiment);
        
        V = [W zfunc(W)];
        ns(i) = size(W,1);
        E = edges(F);
        edgelens(i) = mean(normrow(V(E(:,2),:) - V(E(:,1),:)));
        
        if i==1
            %Print wireframe
            clf;
            set(gcf,'WindowStyle','normal');
            hold on;
            t = {};
            off = 1;
            params = {'FaceColor','flat', 'FaceLighting','flat', ...
                'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
                'AmbientStrength',0.6, 'LineWidth',2};
            t{end+1} = tsurf(F, V, 'CData', zeros(size(V,1),1), params{:});
            l = light('Position',[0 -10 10],'Style','infinite');
            colormap([0.8 0.8 0.8]);
            colorbar;
            switch convexp
                case 'omega'
                    view([-50 40]);
                otherwise
                    view(3);
            end
            camproj('persp');
            axis equal;
            axis off;
            set(gcf,'Color','w');
            set(gca,'Pos',[0 0 1 1]);
            set(gcf,'Pos',[0 0 1920 1000]);
            title([convexperiment ' wireframe']);
            saveas(gcf, ['forwardexperiments-wireframe-' convexperiment '.png']);
        end
        
        n = size(V,1);
        u = f(W);
        Q = curved_hessian(V,F);
        vals(i) = u' * Q * u;
        
        plotV = V;
        plotF = F;
        
        % Refine mesh
        [W,F] = loop(W,F);
    end
    
    %Metric of a monge patch
    switch convexperiment
        case 'spiral'
            zfunc = zfuncnopl;
    end
    syms x1 x2;
    assume(x1,'real');
    assume(x2,'real');
    metric = [1+diff(zfunc([x1 x2]),x1)^2 ...
        diff(zfunc([x1 x2]),x1)*diff(zfunc([x1 x2]),x2); ...
        diff(zfunc([x1 x2]),x2)*diff(zfunc([x1 x2]),x1) ...
        1+diff(zfunc([x1 x2]),x2)^2];
    invmetric = inv(metric);
    df = [diff(f([x1 x2]),x1); diff(f([x1 x2]),x2)];%in coords x1,x2
    dS = sqrt(abs(det(metric))); %the surface element
    %Gaussian curvature of a Monge patch
    kappa = (diff(zfunc([x1 x2]),x1,2)*diff(zfunc([x1 x2]),x2,2) ...
        - diff(diff(zfunc([x1 x2]),x1),x2)^2) / ...
        (1 + diff(zfunc([x1 x2]),x1)^2 + diff(zfunc([x1 x2]),x2)^2)^2;
    curvatureTerm = kappa * df'*invmetric*df;
    %Christoffel symbols of a Monge patch
    Chris(:,:,1) = [diff(diff(zfunc([x1 x2]),x1),x1) ...
        diff(diff(zfunc([x1 x2]),x1),x2); ...
        diff(diff(zfunc([x1 x2]),x2),x1) ...
        diff(diff(zfunc([x1 x2]),x2),x2)] * diff(zfunc([x1 x2]),x1) ...
        /(1+diff(zfunc([x1 x2]),x1)^2+diff(zfunc([x1 x2]),x2)^2);
    Chris(:,:,2) = [diff(diff(zfunc([x1 x2]),x1),x1) ...
        diff(diff(zfunc([x1 x2]),x1),x2); ...
        diff(diff(zfunc([x1 x2]),x2),x1) ...
        diff(diff(zfunc([x1 x2]),x2),x2)] * diff(zfunc([x1 x2]),x2) ...
        /(1+diff(zfunc([x1 x2]),x1)^2+diff(zfunc([x1 x2]),x2)^2);
    nabladf = [diff(df(1),x1) - dot(df,squeeze(Chris(1,1,:))) ...
        diff(df(1),x2) - dot(df,squeeze(Chris(1,2,:))); ...
        diff(df(2),x1) - dot(df,squeeze(Chris(2,1,:))) ...
        diff(df(2),x2) - dot(df,squeeze(Chris(2,2,:)))];
    nabladfdotnabladf = trace(invmetric*nabladf'*invmetric*nabladf);
    Hessintegrand = (nabladfdotnabladf + curvatureTerm) * dS;
    %Hessintegrand = f([x1,x2])^2 * dS;
    %Hessintegrand = nabladfdotnabladf * dS;
    Hessintegrand_fct = matlabFunction(Hessintegrand);
    clear x1 x2;
    
    %Perform quadrature, upsample for accuracy
    M = massmatrix(W,F);
    hess_en_exact = sum(M*bsxfun(Hessintegrand_fct,W(:,1),W(:,2)));
    
    %Compute error
    errs = abs(vals - hess_en_exact);
    %     exps = polyfit(log(ns),log(err),1);
    %     deg_of_c = exps(1);
    
    minel = min(minel, min(edgelens));
    maxel = max(maxel, max(edgelens));
    all_edgelengths{convexp} = edgelens;
    minn = min(minn, min(ns));
    maxn = max(maxn, max(ns));
    all_ns{convexp} = ns;
    all_errs{convexp} = errs;
    
    
    % Plot in terms of number of vertices
    
    close;
    figure;
    set(gcf,'WindowStyle','normal');
    loglog(ns, errs, '-ob', ...
        ns,ns.^(-0.5), '--k');
    title([convexperiment ' convergence plot (verts)']);
    xlabel('number of vertices');
    ylabel('error');
    legend('forward error of fct', ...
        'n^{-0.5}');
    
    saveas(gcf, ['forwardexperiments-' convexperiment '-plot-verts.png']);
    saveas(gcf, ['forwardexperiments-' convexperiment '-plot-verts.eps']);
    
    
    % Plot in terms of average edge length
    
    close;
    figure;
    set(gcf,'WindowStyle','normal');
    loglog(edgelens, errs, '-ob', ...
        edgelens, edgelens, '--k');
    set(gca, 'XDir','reverse');
    title([convexperiment ' convergence plot (edgelens)']);
    xlabel('average edge length');
    ylabel('error');
    legend('forward error of fct', ...
        'h');
    
    saveas(gcf, ['forwardexperiments-' convexperiment '-plot-edgelens.png']);
    saveas(gcf, ['forwardexperiments-' convexperiment '-plot-edgelens.eps']);
    
    
    close;
    figure;
    hold on;
    t = {};
    off = 1;
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
        'AmbientStrength',0.6, 'EdgeColor','none'};
    t{end+1} = tsurf(plotF, plotV, 'CData', u, params{:});
    l = light('Position',[0 -10 10],'Style','infinite');
    %add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
    nc = 50;
    CM = flipud(cbrewer('RdPu'  ,nc));
    colormap(CM);
    %add_isolines(t,'LineWidth',2);
    colorbar;
    switch convexp
        case 'omega'
            view([-50 40]);
        otherwise
            view(3);
    end
    camproj('persp');
    axis equal;
    axis off;
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    %set(gcf,'Pos',[0 0 1920 1000]);
    title([convexperiment ' solution']);
    
    saveas(gcf, ['forwardexperiments-' convexperiment '.png']);
    
end


close;
figure;
loglog(all_ns{1}, all_errs{1}, '-ob', ...
    all_ns{2}, all_errs{2}, '-or', ...
    all_ns{3}, all_errs{3}, '-og', ...
    all_ns{4}, all_errs{4}, '-ob', ...
    all_ns{5}, all_errs{5}, '-om', ...
    all_ns{6}, all_errs{6}, '-oc', ...
    linspace(minn,maxn,10),linspace(minn,maxn,10).^(-0.5), '--k');
title(['all convergence plot (verts)']);
xlabel('number of vertices');
ylabel('error');
legend('blob','omega','arm','spiral','churchdome','octopus', ...
    'n^{-0.5}');
set(gcf,'Pos',[0 0 1920 1000]);

saveas(gcf, ['forwardexperiments-all-plot-verts.png']);
saveas(gcf, ['forwardexperiments-all-plot-verts.eps']);


%Plot error in terms of average edge length

close;
figure;
loglog(all_edgelengths{1}, all_errs{1}, '-ob', ...
    all_edgelengths{2}, all_errs{2}, '-or', ...
    all_edgelengths{3}, all_errs{3}, '-og', ...
    all_edgelengths{4}, all_errs{4}, '-ob', ...
    all_edgelengths{5}, all_errs{5}, '-om', ...
    all_edgelengths{6}, all_errs{6}, '-oc', ...
    linspace(minel,maxel,10),linspace(minel,maxel,10), '--k');
set(gca, 'XDir','reverse');
title(['all convergence plot (edgelens)']);
xlabel('average edge length');
ylabel('error');
legend('blob','omega','arm','spiral','churchdome','octopus', ...
    'h');
set(gcf,'Pos',[0 0 1920 1000]);

saveas(gcf, ['forwardexperiments-all-plot-edgelens.png']);
saveas(gcf, ['forwardexperiments-all-plot-edgelens.eps']);
