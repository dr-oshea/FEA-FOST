function post_processing_NFEA(input, CASE, storeDisps, storeLoads, numLoadSteps, F, PK2, tau, sig, PK1, dofPerNode, nodesPerElem, elementType, Lx, Ly, Lz, meshDivisions)

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

%% plot Stress v Stretch

fprintf('\n   Plotting stress vs stretch curve...');
fprintf('\n');

nGradientComponents = dofPerNode * dofPerNode;
xMin = min(input.ND(:,2));
yMin = min(input.ND(:,3));
spacing = (Lx-xMin)/meshDivisions;

% node to plot stress v stretch at
if elementType == 1
    node = find(input.ND(:,2)==Lx/2 & input.ND(:,3)==Ly/2);
elseif elementType == 2
    % node = find(input.ND(:,2)==Lx/2 & input.ND(:,3)==Ly/2 & input.ND(:,4)==Lz/2);
    node = find(input.ND(:,2)==Lx & input.ND(:,3)==Ly & input.ND(:,4)==Lz);
end

storeF = zeros(numLoadSteps,nGradientComponents);
storeTau = zeros(numLoadSteps,nGradientComponents);
storePK1 = zeros(numLoadSteps,nGradientComponents);
storePK2 = zeros(numLoadSteps,nGradientComponents);
storeSig = zeros(numLoadSteps,nGradientComponents);

connectedElements = zeros(nodesPerElem,2);
for load = 1 : numLoadSteps
    
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

    % Get Cauchy components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = sig(element,nGradientComponents*(localNode-1)+component,load);
        end
        storeSig(load,component) = mean(values);
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
xlabel('\lambda_{1}')
ylabel('\tau_{xx}  [MPa]')

fig = figure;
set(fig,'Name','Stress vs Stretch (YY Component)','NumberTitle','off')
plot(storeF(:,2),storeTau(:,2),'-o');
xlabel('\lambda_{1}')
ylabel('\tau_{yy}  [MPa]')

fig = figure;
hold on
set(fig,'Name','Tensile Stress vs Stretch (2nd Piola Kirchhoff)','NumberTitle','off')
plot(0.5*(storeF(:,1).^2-1),storePK2(:,1),'-o');
plot(0.5*(storeF(:,2).^2-1),storePK2(:,2),'-x');
xlabel('E_{xx}, E_{yy}')
ylabel('S_{xx}, S_{yy} [MPa]')
hold off

fig = figure;
[~,ind] = max(storeF(end,4:end));
set(fig,'Name','Shear Stress vs Stretch (1st Piola Kirchhoff)','NumberTitle','off')
plot(storeF(:,3+ind),storePK1(:,3+ind),'-o');
xlabel('\gamma')
ylabel('PK1 Shear Stress, P [MPa]')

fig = figure;
[~,ind] = max(storeF(end,4:end));
set(fig,'Name','Shear Stress vs Stretch (Cauchy)','NumberTitle','off')
plot(storeF(:,3+ind),storeSig(:,3+ind),'-o');
xlabel('\gamma')
ylabel('Cauchy Shear Stress, \sigma [MPa]')

fig = figure;
set(fig,'Name','Lateral Stretch','NumberTitle','off')
plot(storeF(:,1),storeF(:,2),'-o');
xlabel('\lambda_{1}')
ylabel('\lambda_{2}')

fprintf('\n   Max lambda1 = %.2f',max(storeF(:,1)));
if elementType == 2
    fprintf('\n   Max gamma = %.2f',max(storeF(:,6)));
end


%% Plot stretch vs stress, calculated using force v displacement

% Uniaxial Tension Plotting
if CASE == 2 || CASE == 6 || CASE == 5
    endNodes = find(input.ND(:,2)==Lx);
    endDofX = dofPerNode.*endNodes-(dofPerNode-1).*ones(numel(endNodes),1);

    endStress = zeros(numLoadSteps,1); endStretch = zeros(numLoadSteps,1);
    for load = 1:numLoadSteps
        endForce = sum(storeLoads(endDofX,load));
        endStress(load) = (2*endForce) / (2*Ly*Lz);    % tau = F/A

        endDisp = mean(storeDisps(endDofX,load));
        endStretch(load) = (2*endDisp) / (Lx) + 1;    % lam = Q/L + 1
    end

    fig = figure;
    set(fig,'Name','Calculated End Stress v Stretch');
    plot(endStretch,endStress);
    xlabel('Stretch = 1 + 2*ave(Qx)/xdim')
    ylabel('P11 = 2*sum(Fx)/2*ydim')
end

% Ciarletta-type Shear Experiment
if CASE == 3
    endNodes = find(input.ND(:,3)==Ly);
    endDofX = dofPerNode.*endNodes-(dofPerNode-1).*ones(numel(endNodes),1);

    endStress = zeros(numLoadSteps,1); endStretch = zeros(numLoadSteps,1);
    for load = 1:numLoadSteps
        endForce = sum(storeLoads(endDofX,load));
        endStress(load) = (endForce) / (Lx*Lz);    % tau = F/A

        endDisp = mean(storeDisps(endDofX,load));
        endStretch(load) = (endDisp) / (Ly);    % gam = Q/h
    end

    fig = figure;
    set(fig,'Name','Calculated End Stress v Stretch');
    hold on;
    plot(endStretch,endStress);
    xlabel('gamma = ave(Qx)/ydim')
    ylabel('P12 = sum(Fx)/(xdim*zdim)')
    % xlim([0 0.2])
    % ylim([0 0.045])
end

hold off

%% STRESS COUNTOURS
PLOT = 0;
if elementType == 1 && PLOT ~= 0

xLoc = reshape(input.ND(:,2),[meshDivisions+1,numel(input.ND(:,1))/(meshDivisions+1)]);
yLoc = reshape(input.ND(:,3),[meshDivisions+1,numel(input.ND(:,1))/(meshDivisions+1)]);

% CALCULATE ALL ELEMENT VALUES OF COMPONENTS
elementsToAnalyse = input.EL(:,1);
[F,C,E0,En,S,tau,P,q] = postNonlin_2D(input,storeDisps,storeLoads,elementsToAnalyse);

% Specify which load steps to plot contours at
loadStepsToAnalyse = [1; 2; floor(numLoadSteps/2); (numLoadSteps-10:1:numLoadSteps)'];
screenSize = get(groot,'ScreenSize');

for i = 1:numel(loadStepsToAnalyse)
    loadStep = loadStepsToAnalyse(i);
    lam1 = zeros(nNodes,1);
    lam2 = zeros(nNodes,1);
    tau11 = zeros(nNodes,1);
    tau22 = zeros(nNodes,1);
    tau12 = zeros(nNodes,1);
    E011 = zeros(nNodes,1);
    E022 = zeros(nNodes,1);
    En11 = zeros(nNodes,1);
    En22 = zeros(nNodes,1);
    
    for a = 1:nNodes
        node = input.ND(a,1);
    
        % Get elements connected to specified node
        element1 = find(input.EL(:,1+2) == node);
        element2 = find(input.EL(:,2+2) == node);
        element3 = find(input.EL(:,3+2) == node);
        element4 = find(input.EL(:,4+2) == node);
    
        % Get F11 from each element (lambda1)
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = F(element1,numGradientComponents*(1-1)+1,loadStep); end
        if isempty(element2)==0; value2 = F(element2,numGradientComponents*(2-1)+1,loadStep); end
        if isempty(element3)==0; value3 = F(element3,numGradientComponents*(3-1)+1,loadStep); end
        if isempty(element4)==0; value4 = F(element4,numGradientComponents*(4-1)+1,loadStep); end
        lam1(a) = mean([value1 value2 value3 value4]);
        
        % Get F22 from each element (lambda2)
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = F(element1,numGradientComponents*(1-1)+2,loadStep); end
        if isempty(element2)==0; value2 = F(element2,numGradientComponents*(2-1)+2,loadStep); end
        if isempty(element3)==0; value3 = F(element3,numGradientComponents*(3-1)+2,loadStep); end
        if isempty(element4)==0; value4 = F(element4,numGradientComponents*(4-1)+2,loadStep); end
        lam2(a) = mean([value1 value2 value3 value4]);

        % Get E011 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = E0(element1,numGradientComponents*(1-1)+1,loadStep); end
        if isempty(element2)==0; value2 = E0(element2,numGradientComponents*(2-1)+1,loadStep); end
        if isempty(element3)==0; value3 = E0(element3,numGradientComponents*(3-1)+1,loadStep); end
        if isempty(element4)==0; value4 = E0(element4,numGradientComponents*(4-1)+1,loadStep); end
        E011(a) = mean([value1 value2 value3 value4]);
        
        % Get E022 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = E0(element1,numGradientComponents*(1-1)+2,loadStep); end
        if isempty(element2)==0; value2 = E0(element2,numGradientComponents*(2-1)+2,loadStep); end
        if isempty(element3)==0; value3 = E0(element3,numGradientComponents*(3-1)+2,loadStep); end
        if isempty(element4)==0; value4 = E0(element4,numGradientComponents*(4-1)+2,loadStep); end
        E022(a) = mean([value1 value2 value3 value4]);
        
        % Get En11 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = En(element1,numGradientComponents*(1-1)+1,loadStep); end
        if isempty(element2)==0; value2 = En(element2,numGradientComponents*(2-1)+1,loadStep); end
        if isempty(element3)==0; value3 = En(element3,numGradientComponents*(3-1)+1,loadStep); end
        if isempty(element4)==0; value4 = En(element4,numGradientComponents*(4-1)+1,loadStep); end
        En11(a) = mean([value1 value2 value3 value4]);
        
        % Get En22 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = En(element1,numGradientComponents*(1-1)+2,loadStep); end
        if isempty(element2)==0; value2 = En(element2,numGradientComponents*(2-1)+2,loadStep); end
        if isempty(element3)==0; value3 = En(element3,numGradientComponents*(3-1)+2,loadStep); end
        if isempty(element4)==0; value4 = En(element4,numGradientComponents*(4-1)+2,loadStep); end
        En22(a) = mean([value1 value2 value3 value4]);
        
        % Get tau11 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = tau(element1,numGradientComponents*(1-1)+1,loadStep); end
        if isempty(element2)==0; value2 = tau(element2,numGradientComponents*(2-1)+1,loadStep); end
        if isempty(element3)==0; value3 = tau(element3,numGradientComponents*(3-1)+1,loadStep); end
        if isempty(element4)==0; value4 = tau(element4,numGradientComponents*(4-1)+1,loadStep); end
        tau11(a) = mean([value1 value2 value3 value4]);
        
        % Get tau22 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = tau(element1,numGradientComponents*(1-1)+2,loadStep); end
        if isempty(element2)==0; value2 = tau(element2,numGradientComponents*(2-1)+2,loadStep); end
        if isempty(element3)==0; value3 = tau(element3,numGradientComponents*(3-1)+2,loadStep); end
        if isempty(element4)==0; value4 = tau(element4,numGradientComponents*(4-1)+2,loadStep); end
        tau22(a) = mean([value1 value2 value3 value4]);
        
        % Get tau12 from each element
        value1=[]; value2=[]; value3=[]; value4=[];
        if isempty(element1)==0; value1 = tau(element1,numGradientComponents*(1-1)+3,loadStep); end
        if isempty(element2)==0; value2 = tau(element2,numGradientComponents*(2-1)+3,loadStep); end
        if isempty(element3)==0; value3 = tau(element3,numGradientComponents*(3-1)+3,loadStep); end
        if isempty(element4)==0; value4 = tau(element4,numGradientComponents*(4-1)+3,loadStep); end
        tau12(a) = mean([value1 value2 value3 value4]);
        
    end
    
    fig = figure('OuterPosition',[1 1 screenSize(3) screenSize(4)]);
    str = sprintf('Contour Plots: Load Step #%i',loadStep);
    set(fig,'Name',str,'NumberTitle','off');

    subplot(2,4,1)
    Z = reshape(lam1,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('Lambda 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(lam1)));
    xlabel(str)
    colorbar
    
    subplot(2,4,2)
    Z = reshape(lam2,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('Lambda 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(lam2)));
    xlabel(str)
    colorbar
    
    subplot(2,4,3)
    Z = reshape(tau11,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('tau 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(tau11)));
    xlabel(str)
    colorbar
    
    subplot(2,4,4)
    Z = reshape(tau22,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('tau 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(tau22)));
    xlabel(str)
    colorbar
    
    subplot(2,4,5)
    Z = reshape(E011,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('E0 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(E011)));
    xlabel(str)
    colorbar
    
    subplot(2,4,6)
    Z = reshape(E022,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('E0 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(E022)));
    xlabel(str)
    colorbar
    
    subplot(2,4,7)
    Z = reshape(En11,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('En 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(En11)));
    xlabel(str)
    colorbar
    
    subplot(2,4,8)
    Z = reshape(En22,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('En 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(En22)));
    xlabel(str)
    colorbar
    
    figure
    Z = reshape(tau12,[meshDivisions+1,meshDivisions+1]);
    contourf(xLoc,yLoc,Z,'ShowText','on');
    title('tau 12')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(tau12)));
    xlabel(str)
    colorbar
    
    
end
end

end