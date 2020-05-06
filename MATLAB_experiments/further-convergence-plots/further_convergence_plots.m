addpath ../../../cpp/build/applications/Mex
addpath ..

convexperiments = {'plane','cylinder','lilium','bob','hand','cheburashka'};
%convexperiments = {'hand'};


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
        case 'plane'
            nsubdivs = 9;
            [V,F] = create_regular_grid(4, 4);
            V = [(V - [0.5 0.5] ) * 2 zeros(size(V,1),1)];
            phix = 0;
        case 'cylinder'
            nsubdivs = 9;
            [X,Y,Z] = cylinder(ones(4,1),6);
            [Q,V] = surf2patch(X,Y,Z);
            F = [Q(:,1) Q(:,2) Q(:,3); Q(:,1) Q(:,3) Q(:,4)];
            [V,~,IM] = remove_duplicate_vertices(V,eps);
            F = IM(F);
            cylproj = @(x) [(1.5 - x(:,3).^2) .* x(:,1:2) ./ ...
                normrow(x(:,1:2)) x(:,3)];
            V = cylproj(V);
            phix = 0;
        case 'lilium'
            nsubdivs = 8;
            [V,F] = readOBJ('lilium.obj');
            phix = 0;
        case 'bob'
            nsubdivs = 8;
            [V,F] = readOBJ('bob.obj');
            phix = pi/2;
        case 'hand'
            nsubdivs = 7;
            [V,F] = readOBJ('hand.obj');
            phix = 0;
        case 'cheburashka'
            nsubdivs = 7;
            [V,F] = readOBJ('cheburashka.obj');
            phix = 0;
    end
    
    R = [1 0 0; 0 cos(phix) -sin(phix); 0 sin(phix) cos(phix)];
    V = V*R';
    
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
    if(strcmp(convexperiment,'hand'))
        view([18, 20]);
    else
        view(3);
    end
    camproj('persp');
    axis equal;
    axis off;
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    title([convexperiment ' wireframe']);
    saveas(gcf, ['further-experiments-wireframe-' convexperiment '.png']);
    
    S = {speye(size(V,1), size(V,1))};
    zs = {};
    ns = zeros(1,nsubdivs);
    edgelens = zeros(1,nsubdivs);
    
    for i=1:nsubdivs
        fprintf('computing upsample %d for experiment %s\n', i, ...
            convexperiment);
        
        switch convexperiment
            case 'plane'
                bl = find(V(:,1)<-1+1e-12);
                br = find(V(:,1)>1-1e-12);
                bc = [ones(size(bl,1),1); zeros(size(br,1),1)];
                b = [bl; br];
            case 'cylinder'
                bb = find(V(:,3)<1e-12);
                bt = find(V(:,3)>1-1e-12);
                bc = [ones(size(bb,1),1); zeros(size(bt,1),1)];
                b = [bb; bt];
            case 'lilium'
                b = unique(outline(F));
                bc = cos(2*atan2(V(b,2), V(b,1)));
            case 'bob'
                b = unique(outline(F));
                bl = b(V(b,2)<0);
                br = b(V(b,2)>0);
                b = [bl; br];
                bc = [ones(size(bl)); zeros(size(br))];
            case 'hand'
                b = unique(outline(F));
                bc = -(V(b,1).^2 + 0.1*V(b,2) + sin(pi*V(b,3)));
            case 'cheburashka'
                b = unique(outline(F));
                bc = sin(2*atan2(V(b,2), V(b,1)));
        end
        
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
            
            switch convexperiment
                case 'plane'
                    % Set linf-norm of boundary to 1.
                    newb = unique(outline(F));
                    V(newb,:) = V(newb,:) ./ max(abs(V(newb,:)),[],2);
                case 'cylinder'
                    V = cylproj(V);
                    newb = unique(outline(F));
                    V(newb,3) = V(newb,3)>0.5;
                case 'lilium'
                    newb = unique(outline(F));
                    V(newb,3) = 0;
                case 'bob'
                    b = unique(outline(F));
                    V(b,1) = 0;
                case 'hand'
                    b = unique(outline(F));
                    V(b,3) = 0;
                case 'cheburashka'
                    b = unique(outline(F));
                    V(b,3) = 0;
            end
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
    
    exps = polyfit(log(ns),log(err),1);
    deg_of_c = exps(1);
    
    minel = min(minel, min(edgelens));
    maxel = max(maxel, max(edgelens));
    all_edgelengths{convexp} = edgelens;
    minn = min(minn, min(ns));
    maxn = max(maxn, max(ns));
    all_ns{convexp} = ns;
    all_errs{convexp} = err;
    
    %Plot error in terms of number of vertices
    
        close;
        figure;
        set(gcf,'WindowStyle','normal');
        loglog(ns, err, '-ob', ...
            ns,ns.^(-0.5), '--k');
        title([convexperiment ' convergence plot (verts)']);
        xlabel('number of vertices');
        ylabel('error');
        legend('l2 error of interpolated data', ...
            'n^{-0.5}');
    
        saveas(gcf, ['further-experiments-' convexperiment '-plot-verts.png']);
        saveas(gcf, ['further-experiments-' convexperiment '-plot-verts.eps']);
    
    
        %Plot error in terms of average edge length
    
        close;
        figure;
        set(gcf,'WindowStyle','normal');
        loglog(edgelens, err, '-ob', ...
            edgelens,edgelens, '--k');
        set(gca, 'XDir','reverse');
        title([convexperiment ' convergence plot (edgelens)']);
        xlabel('average edge length');
        ylabel('error');
        legend('l2 error of interpolated data', ...
            'h');
    
        saveas(gcf, ['further-experiments-' convexperiment '-plot-edgelens.png']);
        saveas(gcf, ['further-experiments-' convexperiment '-plot-edgelens.eps']);
    
    
    clf;
    set(gcf,'WindowStyle','normal');
    hold on;
    t = {};
    off = 1;
    params = {'FaceColor','interp', 'FaceLighting','gouraud', ...
        'SpecularStrength',0.1, 'DiffuseStrength',0.8, ...
        'AmbientStrength',0.6, 'EdgeColor','none'};
    t{end+1} = tsurf(F, V, 'CData', zs{end}, params{:});
    l = light('Position',[0 -10 10],'Style','infinite');
    %add_shadow(t,l,'Color',[0.7 0.7 0.7],'Fade','infinite');
    nc = 50;
    CM = flipud(cbrewer('GnBu'  ,nc));
    colormap(CM);
    %add_isolines(t,'LineWidth',2);
    colorbar;
    if(strcmp(convexperiment,'hand'))
        view([18, 20]);
    else
        view(3);
    end
    camproj('persp');
    axis equal;
    axis off;
    set(gcf,'Color','w');
    set(gca,'Pos',[0 0 1 1]);
    set(gcf,'Pos',[0 0 1920 1000]);
    title([convexperiment ' solution']);
    
    saveas(gcf, ['further-experiments-' convexperiment '.png']);
    
end



close;
figure;
loglog(all_ns{1}, all_errs{1}, '-ob', ...
    all_ns{2}, all_errs{2}, '-or', ...
    all_ns{3}, all_errs{3}, '-og', ...
    all_ns{4}, all_errs{4}, '-oy', ...
    all_ns{5}, all_errs{5}, '-om', ...
    all_ns{6}, all_errs{6}, '-oc', ...
    linspace(minn,maxn,10),linspace(minn,maxn,10).^(-0.5), '--k');
title(['all convergence plot (verts)']);
xlabel('number of vertices');
ylabel('error');
legend('plane','cylinder','lilium','bob','hand','cheburashka', ...
    'n^{-0.5}');
set(gcf,'Pos',[0 0 1920 1000]);

saveas(gcf, ['further-experiments-all-plot-verts.png']);
saveas(gcf, ['further-experiments-all-plot-verts.eps']);


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
legend('plane','cylinder','lilium','bob','hand','cheburashka', ...
    'h');
set(gcf,'Pos',[0 0 1920 1000]);

saveas(gcf, ['further-experiments-all-plot-edgelens.png']);
saveas(gcf, ['further-experiments-all-plot-edgelens.eps']);