%% Function to do post-processing plots - Kiru

s = 'saved_location.txt';
dPath = which (s);
cd (dPath(1:end-length(s)));

load("input.mat", "input")
load("displacements.mat", "storeDisps")
load("loads.mat", "storeLoads")

dofPerNode = 3;
nodesPerElem = 8;
elementType = 4;
numLoadSteps = size(storeDisps,2);

%%                      Post-processing of Results

% This section runs an algorithm which returns arrays storing various
% deformation and stress tensors at each node (plus centre) of each
% element. This can later be used for plotting results

fprintf('\n   Commencing Post Processing...')

% Choose elements to analyse
elementsToAnalyse = input.EL(:,1);    % all elements

% Run post processing (returns deformation and stress vectors)
if elementType == 1
    [F,~,~,~,PK2,tau,PK1,~] = postNonlin_2D(input,storeDisps,storeLoads,elementsToAnalyse);
elseif elementType == 2 || elementType == 4
    [F,~,~,~,PK2,tau,sig,PK1,~] = postNonlin_3D(input,storeDisps,storeLoads,elementsToAnalyse);
end

%%                      Plotting of Results

% End game - would be nice to plot changing displacement/stress contour
% plots as load step progresses (use app designer)

fprintf('\n   Plotting deformed configurations at each load step...');
fprintf('\n');

coords0 = input.ND;
def = [zeros(size(input.ND,1),1) reshape(storeDisps(:,end),dofPerNode,size(input.ND,1))'] + coords0;
xMin = min([min(input.ND(:,2)),min(def(:,2))]);
xmax = max([max(input.ND(:,2)),max(def(:,2))]);
yMin = min([min(input.ND(:,3)),min(def(:,3))]);
ymax = max([max(input.ND(:,3)),max(def(:,3))]);
if elementType == 2
    zmin = min([min(input.ND(:,4)),min(def(:,4))]);
    zmax = max([max(input.ND(:,4)),max(def(:,4))]);
end

M(numLoadSteps) = struct('cdata',[],'colormap',[]);

% Plot Deformed Shapes
for load = numLoadSteps:-max(factor(numLoadSteps)):1

    Q = storeDisps(:,load);

    Qmap = reshape(Q,dofPerNode,size(input.ND,1))';

    input.ND(:,2:end) = coords0(:,2:end) + Qmap;

    fig = figure;
    fig.Visible = 'off';
    str = sprintf('FEM Solution: Load Step #%i',load);
    set(fig,'Name',str,'NumberTitle','off');
    plotDiscretisation(input,0,elementType);
    xlim([xMin xmax])
    ylim([yMin ymax])
    if elementType == 2
        zlim([zmin zmax]);
    end
    drawnow
    M(load) = getframe;

end
input.ND = coords0; % reset the input data to original coordinates
fig.Visible = 'on';

% Produce animation of deformed shape
movie(M(numLoadSteps:-max(factor(numLoadSteps)):1),10)

%% plot contour shapes
load_step = numLoadSteps;% uses last load step.
[nodal_F, nodal_S] = convert_element_stress_to_nodal_stress(F, sig, input, load_step);

nonlinear_postprocessing_plotting(input, storeDisps(:,end), nodal_S, 3)

%% plot Stress v Stretch

fprintf('\n   Plotting stress vs stretch curve...');
fprintf('\n');

nGradientComponents = dofPerNode * dofPerNode;

% node to plot stress v stretch at
if elementType == 1
    node = find(input.ND(:,2)==Lx/2 & input.ND(:,3)==Ly/2);
elseif elementType == 2
    % node = find(input.ND(:,2)==Lx/2 & input.ND(:,3)==Ly/2 & input.ND(:,4)==Lz/2);
    node = find(input.ND(:,2)==Lx & input.ND(:,3)==Ly & input.ND(:,4)==Lz);
elseif elementType == 4
    % get node on top circle, outer-surface.
    Xloc = 0;
    Yloc = max(input.ND(:,3));
    Zloc = 0;
    tol = 10^-2;
    node = find(abs(input.ND(:,2)-Xloc) <= tol & abs(input.ND(:,3)-Yloc) <= tol & abs(input.ND(:,4)-Zloc) <= tol);
    node = 97;   %USE THIS TO SELECT NODE
end



storeF = zeros(numLoadSteps,nGradientComponents);
storeTau = zeros(numLoadSteps,nGradientComponents);
storePK1 = zeros(numLoadSteps,nGradientComponents);
storePK2 = zeros(numLoadSteps,nGradientComponents);

connectedElements = zeros(nodesPerElem,2);
for load = 1 : numLoadSteps
    

    % change coordinates from XYZ to FSN
    convert_to_FSN = 1;
    element_F = F(:,:,load);
    element_S = sig(:,:,load);
    if convert_to_FSN == 1
        element_F_new = element_F*0;
        element_S_new = element_S*0;
        for elem = 1:size(input.EL,1)
            f = input.FIBRES.f(elem,:);
            s = input.FIBRES.s(elem,:);
            n = input.FIBRES.n(elem,:);
            Q = [f; s; n]';
            for nodeX = 1:nodesPerElem
                Fvec = element_F(elem,(nodeX-1)*nGradientComponents+1:nodeX*nGradientComponents);
                Fmat = V2M(Fvec');
                F_fsn = Q.'*Fmat*Q;
    
                S_fsn = Get_S_CC(F_fsn,input.model,input.FIBRES);
                sig_fsn = F_fsn*S_fsn*F_fsn.';
                element_F_new(elem,(nodeX-1)*nGradientComponents+1:nodeX*nGradientComponents) = M2V(F_fsn)';
                element_S_new(elem,(nodeX-1)*nGradientComponents+1:nodeX*nGradientComponents) = M2V(sig_fsn)';
            end
        end
        F(:,:,load) = element_F_new;
        tau(:,:,load) = element_S_new;
    end


    % This could be improved by using another loop and arrays..
    i = 0;
    for localNode = 1 : nodesPerElem
        iElement = find(input.EL(:,2+localNode) == node);
        if isempty(iElement) ~= 1
            i = i + 1;
            connectedElements(i,:) = [iElement, localNode];
        end
    end
    numConnected = i;

    % Get F components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = F(element,nGradientComponents*(localNode-1)+component,load);
        end
        storeF(load,component) = mean(values);
    end

    % Get tau components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = tau(element,nGradientComponents*(localNode-1)+component,load);
        end
        storeTau(load,component) = mean(values);
    end

    % Get PK1 components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = PK1(element,nGradientComponents*(localNode-1)+component,load);
        end
        storePK1(load,component) = mean(values);
    end

    % Get PK2 components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = PK2(element,nGradientComponents*(localNode-1)+component,load);
        end
        storePK2(load,component) = mean(values);
    end

           
end    
    
fig = figure;
set(fig,'Name','Stress vs Stretch (XX Component)','NumberTitle','off')
plot(storeF(:,1),storeTau(:,1),'-o');
xlabel('\lambda_{F}')
ylabel('\sigma_{F}  [MPa]')

fig = figure;
set(fig,'Name','Stress vs Stretch (YY Component)','NumberTitle','off')
plot(storeF(:,2),storeTau(:,2),'-o');
xlabel('\lambda_{S}')
ylabel('\sigma_{SS}  [MPa]')

fprintf('\n   Max lambdaF = %.2f',max(storeF(:,1)));
% fprintf('\n   Max gamma = %.2f',max(storeF(:,6)));

