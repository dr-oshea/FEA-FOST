%% FEA_nonlin2D

% This script runs the nonlinear finite element program where the strain
% energy function and corresponding constitutive law is described by the 
% work submitted by A/Prof. Mario Attard and Daniel O'Shea in 2018. The 
% program uses the Initial or Secant Gradient Method to iteratively solve a
% nonlinear hyperelasticity problem in 2D/3D subject to given applied nodal 
% loads and constraints

% Author:   Daniel O'Shea   < d.oshea@unsw.edu.au >
% Created:  16 March 2018

% Update: 20 April 2023
% Modified to include a closed form Tangent stiffness method, following
% publication of Myocardium paper in 2022.

%%                             Preliminary

% Clear the workspace
clc
clear
close all

% Select Case to Study:
% 1 = Pinned bottom edge, tension top edge Y- direction
% 2 = Uniaxial tension in X-Direction
% 3 = Ciarletta shear test
% 4 = Simple shear test
% 5 = Equibiaxial tension
% 6 = Uniaxial Tension 1 Element
% 7 = Pure Shear (Biaxial Testing)
CASE = 3;   % choose boundary conditions from previously set-up cases.
MAT = 4;    % 1 = specify isotropic, 2 = specify orthotropic, 3 = known sets of properties
MIXED = 0;  % use incompressible formulation

fprintf('\n   Commencing Nonlinear FE Analysis Program');
fprintf('\n');

% Start Timer:
tic

%%                           Define the problem

% This section defines the geometry of the continuum, the discretisation,
% the material properties, and the boundary conditions (loads and
% constraints)

% Dimension of problem: Type of Element used
elementType = 3;           % Type of element: 1 = 2D quad, 2 = 3D Brick, 3 = 2d annulus.

% Geometry
Lx = 2*pi;    % [mm]    % maximum angle (should be 2*Pi)
Ly = 4;       % [mm]    % thickness of annulus
Lz = 4;       % [mm]    % inner radius
thetaLayer1 = 0;     % [deg]
thetaLayer2 = 0;     % [deg]
numLayers = 1;        % theta1 and theta2 give orientation of material for case of multiple layers (layers > 1)

% Discretisation
meshDivisions = 5;                     % element divisions along theta
eleAspectRatio = 5/meshDivisions;       % change numerator for element divisions along radius

% Isotropic Material Properties
if MAT == 1
    E = 1.434;
    v = 0.495;
    G = E/(2*(1+v));
    materialProps = [E E E v v v G G G];
end

% Orthotropic Material Properties
if MAT == 2 
    E11 = 0.788;
    E22 = 0.2;
    E33 = 0.2;
    v23 = 0.428;
    v13 = 0.431;
    v12 = 0.431;
    G23 = 0.01;
    G13 = 0.49;
    G12 = 0.49;
    materialProps = [E11 E22 E33 v23 v13 v12 G23 G13 G12];
end

% Material Property Sets
if MAT == 3
    CiarlettaA = [0.488, 0.383, 0.383, 0.428, 0.431, 0.431, 0.15, 0.12, 0.12];
    CiarlettaB = [132.147, 4.021, 4.021, 0.35, 0.268, 0.268, 1.0, 1.3, 1.3];
    HeullietDegautier =  [1.43473346911424;1.43473346911424;1.43473346911424;0.4;0.4;0.4;0.490032668844590;0.490032668844590;0.490032668844590];
    Itskov = [0.4258, 0.0803, 0.0803, 0.249, 0.260, 0.260, 1, 1, 1];
    Treloar = [];
%     Ciarletta2009 = [0.0237;0.0918;0.0059;1.3815;1.1549;0.0016;0.0002;0.0004;0.0004];
    Ciarletta2009 = get_MatProps([[1;1;1].*10^0;[0.01816;0.0746;0.00803;0.0005;0.0007;0.0007]]);
    materialProps = CiarlettaB;
    
end

% Always use MAT = 4 (previous code now redundant)
% Select model, then alter the material constants as desired.
if MAT == 4
    constitutiveModel = 1;

    if constitutiveModel == 1   % Myocardium model (2022)
        L00 = 1000;  % this needs to be orders of magnitude larger than other values (enforces incompressibility)
        LFF = 717;
        LSS = 215;
        LNN = 95;
        M00 = 3;
        MFS = 0;
        kappaF = 0.1738;
        kappaS = 0.1179;
        k = 2;
        input.model = [constitutiveModel, L00, LFF, LSS, LNN, M00, MFS, kappaF, kappaS, k];

    elseif constitutiveModel == 2      % Gultekin model (2017)
        a = 0.4;
        b = 6.55; 
        af = 3.05; 
        bf = 29.05;
        as = 1.25;
        bs = 36.65;
        afs = 0.15;
        bfs = 6.28;
        kappaF = 0.08; 
        kappaS = 0.09;
        input.model = [constitutiveModel, a, b, af, bf, as, bs, afs, bfs, kappaF, kappaS];

    end
end
PLANAR = 2; % 1 = Plane stress, 2 = Plane strain;

% Specified final loading (Force OR Stress used, depending on problem)
numLoadSteps = 50;        % Load steps (increasing takes longer, though more stable solution process)   CHANGED LOAD STEP FROM 250 to 10 to just get it running
finalForce = 0.0;      % [Newtons]
finalStress = 5;      % [should match units of material parameters, kPa, MPa, etc..]

% Solution procedure options
GRAD = 3;       % gradient method used (3 = tangent stiffness)
BROYDEN = 0;    % enable Broyden's quasi-Newton method (not relevant)
PROCMOD = 1;    % enable process modification (not relevant)
QUASI_ITER = 10; % number of iterations before new tangent stiffness calculation
TOL_equilibrium = 10^-3; % tolerance on equilibrium 
TOL_procmod = 10^-2;     % tolerance on process modifications (not relevant)
maxMod = 10;    % maximum number of process "modifications" before exiting
maxIter = 1000;  % maximum number of iterations before exiting 
PLOT = 0;       % choose whether to suppress plotting (not relevant)
PROG = 0;       % display progress bars (this works! - set 0 to not show)

SolverOptions.lineSearchOptions.use = 1;
SolverOptions.lineSearchOptions.method = "Quadratic Expansion";
SolverOptions.lineSearchOptions.tolerance = 0.8;    % set between 0.5 and 0.8
SolverOptions.lineSearchOptions.maxSearch = 15;
SolverOptions.MIXED = 0;

%%                           Generate input arrays
fprintf('\n   Determining Input Information...');

% Boundary conditions for previously set-up cases (can add more)
if CASE == 1
    CCASE = 1;
    FCASE = 1;
elseif CASE == 2
    CCASE = 5; 
    FCASE = 2;
    fprintf('\n   "Uniaxial Tension: X-Direction"');
elseif CASE == 3
    CCASE = 3;
    FCASE = 3;
    fprintf('\n   "Ciarletta Test Problem"');
elseif CASE == 4
    CCASE = 4;
    FCASE = 4;
    fprintf('\n   "Simple Shear Test"');
elseif CASE == 5
    CCASE = 7;
    FCASE = 5;
    fprintf('\n   "Equibiaxial Tension Test"');
elseif CASE == 6
    CCASE = 6;
    FCASE = 6;
    fprintf('\n   "Uniaxial Tension Test: ONE ELEMENT"'); 
elseif CASE == 7
    CCASE = 7;
    FCASE = 2;
    fprintf('\n   "Pure Shear (Biaxial Testing)"');
elseif CASE == 91
    CCASE = 91;    % need to create this case inside input_constraint
    FCASE = 91;    % need to create this case inside input_load 
    fprintf('\n   "Simple Shear (FS mode)"');
elseif CASE == 92
    CCASE = 92;
    FCASE = 92;
    fprintf('\n   "Simple Shear (FN mode)"');
elseif CASE == 93
    CCASE = 93;
    FCASE = 93;
    fprintf('\n   "Simple Shear (SF mode)"');
elseif CASE == 94
    CCASE = 94;
    FCASE = 94;
    fprintf('\n   "Simple Shear (SN mode)"');
elseif CASE == 95
    CCASE = 95;
    FCASE = 95;
    fprintf('\n   "Simple Shear (NF mode)"');
elseif CASE == 96
    CCASE = 96;
    FCASE = 96;
    fprintf('\n   "Simple Shear (NS mode)"');
end
fprintf('\n');

% Set up "input" object: contains info on Nodes, Elements, Loads, and
% Constraints

% Array storing all nodal values
input.ND = input_nodal(meshDivisions,eleAspectRatio,Lx,Ly,Lz,elementType);

% Array storing all element connectivity
input.EL = input_element(input.ND,meshDivisions,eleAspectRatio,numLayers,Lx,Ly,Lz,thetaLayer1,thetaLayer2,elementType);

% Array storing all material properties and material tensors
% MATERIAL PROPERTIES ARE HARD CODED IN in FILE Get_S_CC.m
% [input.C0,input.L0,input.M0] = input_material(MatProps,PLANAR,type);

% Array storing all boundary conditions
input.CON  = input_constraint(input.ND,Lx,Ly,Lz,elementType,CCASE);
input.LOAD = input_load(input.ND,input.EL,Lx,Ly,Lz,meshDivisions,elementType,finalForce,finalStress,FCASE);

% Array storing element information
[dofPerNode, nodesPerElem, strainComps] = input_info(elementType);

% Additional info
nNodes = size(input.ND,1);      % total nodes in defined problem
totalDOFs = nNodes*dofPerNode;  % total degrees of freedom in defined problem
nElems = size(input.EL,1);      % total elements in defined problem


% This is the end of pre-processing.

%% Display details of defined elasticity problem
fprintf('\n   Problem Generated.');
fprintf('\n');
fprintf('\n   Summary of Problem:');
fprintf('\n      Body width  = %.2f mm',Lx);
fprintf('\n      Body height = %.2f mm',Ly);
fprintf('\n      Constriants:');
if CCASE == 1
    fprintf('\n         Bottom edge pinned.');
elseif CCASE == 2
    fprintf('\n         Centre node pinned.');
end
fprintf('\n      Loading:');
if FCASE == 1
    fprintf('\n         Top edge uniform load = %.2f N/mm2',finalStress);
elseif FCASE == 2
    fprintf('\n         Top edge uniform load    =  %.2f N/mm2',finalStress);
    fprintf('\n         Bottom edge uniform load = -%.2f N/mm2',finalStress);
end
fprintf('\n');
fprintf('\n   Summary of Discretisation:');
if elementType == 1
    fprintf('\n      2D Bilinear Quadrilateral Elements');
elseif elementType == 2
    fprintf('\n      3D Hexahedral Elements');
end
fprintf('\n      Number of Nodes:          %i',size(input.ND,1));
fprintf('\n      Number of Elements:       %i',size(input.EL,1));
fprintf('\n      Total Degrees of Freedom: %i',totalDOFs);
fprintf('\n');
fprintf('\n');


%%                           Plot of Discretisation

% This section produces a plot of the discretised mesh and applied boundary
% condtitions (appropriated from open-source code)

fig = figure;
set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
plotDiscretisation(input,1,elementType)

%%                         Determine Initial Gradient

% This section determines the initial gradient of the elasticity problem.
% This is equivalent to the tangent stiffness matrix of each element 
% evaluated in the region of infinitesimal displacement gradients

fprintf('\n   Determining Initial Stiffness Matrix...');

K0 = Assemble_Initial_Stiffness(input,elementType,nodesPerElem,dofPerNode);

fprintf('\n   Completed.');
fprintf('\n');
fprintf('\n   Summary of Initial Stiffness Matrix:');
fprintf('\n      Total size:             %i x %i',dofPerNode*size(input.ND,1),dofPerNode*size(input.ND,1));
fprintf('\n      Total entries:          %i',(dofPerNode*size(input.ND,1))^2);
fprintf('\n      Total non-zero entries: %i',nnz(K0));
fprintf('\n');


%%                   Define Incremental Load Vector

% This section generates the final load vector based on input conditions,
% before using this to determine the global incremental load vector

fprintf('\n   Constructing incremental load vector...');

% Number of load steps
fprintf('\n      Total load steps:     %i',numLoadSteps);

% Final Load Vector
finalLoad = Final_Load_Vector(input);
fprintf('\n      Final load vector constructed');
fprintf('\n      Maximum Load:         %.3f N',max(max(input.LOAD(:,2:end))));

% Initial Load Vector
initialLoad = 0.*finalLoad;

% Incremental load vector
DeltaP = (1/numLoadSteps).*(finalLoad - initialLoad);
fprintf('\n      Incremental load vector determined');
fprintf('\n');
maxInc = max(max(input.LOAD(:,2:end)))/numLoadSteps;

%%                     Commence Iterative Process

% Determine vector of unconstrained degrees of freedom
freeDOFs = Apply_Constraints(input);


%%                         Solve elasticity problem

% This section iteratively estimates the displacement solution at each
% defined load step until specified convergence limits are met

fprintf('\n   SOLVING THE ELASTICITY PROBLEM')
fprintf('\n')
if GRAD == 0
    fprintf('\n   "Initial Stiffness Method" \n');
elseif GRAD == 1
    fprintf('\n   "Secant Stiffness Method"  \n');
elseif GRAD == 2
    fprintf('\n   "Finite Difference Tangent Stiffness Method" \n');    
elseif GRAD == 3
    % Use of tangent stiffness method requires closed form representation
    % of PK2 stress tensor and Tangent elasticity tensor in script
    % "Get_S_CC.m"
    fprintf('\n   "Tangent Stiffness Method" \n');    
end

% Allocate initial displacement and load vectors (zero vectors)
storeDisps = zeros(totalDOFs,numLoadSteps);
Q0   = zeros(totalDOFs,1);
storeLoads = zeros(totalDOFs,numLoadSteps);
Pest = zeros(totalDOFs,numLoadSteps);

% Solution tolerance
fprintf('\n   Solution Tolerance: %.0d',TOL_equilibrium)
fprintf('\n')

% Iterative Process:
maxTar = maxInc;

fprintf('\n   Process Modification:');
fprintf('\n      Eigenvalue Tolerance:  %.0d',TOL_procmod);
fprintf('\n      Maximum modifications: %i',maxMod);
fprintf('\n');

Beta = zeros(maxMod,1);
Phi = zeros(totalDOFs,maxMod);
NEW_TANG = 0;
FAIL = 0;
Kp = zeros(totalDOFs,totalDOFs);

% For each load step
n = 1; newAttempt = 0;
targetLoad = initialLoad;
maxNewAttempts = 3;

while norm(targetLoad) <= norm(finalLoad)
% while n <= numLoadSteps
    fprintf('\n\n   COMMENCING LOAD STEP #%i',n);
    
    % Initial Guess = previous known solution
    if n == 1
        Q = Q0;
        P = initialLoad;
    else
        Q = storeDisps(:,n-1);     % previous known vector
        P = storeLoads(:,n-1)';    % previous known vector
    end
    
    % Targeted Load Vector = previous known + deltaP
    targetLoad = P + DeltaP;
    fprintf('\n      Norm of current target load:   %.5f', norm(targetLoad));
    fprintf('\n      Norm of final target load:     %.5f', norm(finalLoad));

    % Initial Stiffness
    if GRAD == 0 || (GRAD == 2 && n == 1) || (GRAD == 3 && n == 1)
        K = K0;
    
    % Secant Stiffness
    elseif GRAD == 1
        fprintf('\n      Assembling secant stiffness matrix...')
        K = Assemble_Secant_Stiffness(Q,input,elementType,nodesPerElem,dofPerNode);
        
    % Finite Difference Tangent Stiffness
    elseif GRAD == 2 && NEW_TANG == 0
        fprintf('\n      Assembling approximated tangent stiffness matrix...')
        hwb    = waitbar(0,'Assembly process ...');
        set(hwb,'Name','Approximating Tangent Stiffness Matrix');
        indj = 1;
        eta = 10^-3;
        for j = 1:totalDOFs
            if j>(indj*10)
                waitbar(j/totalDOFs,hwb,[num2str(floor(100*j/totalDOFs)) ' % assembly ...']);
                indj=indj+1; 
            end
            if j==totalDOFs
                waitbar(j/totalDOFs,hwb,' assembly finished ');
                disp(['     K assembly: ' num2str(i) ' dofs - complete']);  
            end
            signQ = 1;
            if Q(j) < 0; signQ = -1; end
            h = sqrt(eta)*max(abs(Q(j)),1)*signQ;
            Qj = Q(j);
            Q(j) = Q(j) + h;
            h = Q(j) - Qj;
            Pj = SolveP(Q,input,nodesPerElem);
            for i = 1:totalDOFs
                K(i,j) = (Pj(i)-P(i))/h;
            end
            Q(j) = Qj;
        end
        NEW_TANG = 1;
        close(hwb);
        if n == 30
            hello = 1;
        end
        
    elseif (GRAD == 3 && n > 1)
        fprintf('\n      Assembling TANGENT stiffness matrix...')
        K = Assemble_Tangent_Stiffness(Q,input,elementType,nodesPerElem,dofPerNode,PROG);

    end

    % Output information regarding stiffness matrix in use:
    fprintf('\n      Information regarding stiffness matrix:');
    fprintf('\n         Determinant   = %.2e',det(K(freeDOFs,freeDOFs)));
    fprintf('\n         Condition No. = %.2e',condest(K(freeDOFs,freeDOFs)));
    fprintf('\n         Symmetry      = %.2e',norm(full(K-K.')));
    fprintf('\n');
    
    Kp = K; % Store previous stiffess matrix
   
    % Perform nonlinear iterative process to reach a solution...
    % PROG controls the solution procedure used
    if MIXED == 1
       lambda = 0;
       [Q,P,FAIL,iter,normR] = NonlinearIteration(input,Q,P,targetLoad,K,freeDOFs,totalDOFs,nodesPerElem,TOL_equilibrium,maxIter,maxMod,PROCMOD,QUASI_ITER,FAIL,PROG,SolverOptions);
    else
       [Q,P,FAIL,iter,normR] = NonlinearIteration(input,Q,P,targetLoad,K,freeDOFs,totalDOFs,nodesPerElem,TOL_equilibrium,maxIter,maxMod,PROCMOD,QUASI_ITER,FAIL,PROG,SolverOptions);  
    end
   
    
    % If program completed current load step
    if FAIL == 0

        % Store converged displacement solution
        storeDisps(:,n) = Q;

        % Store target load vector
        storeLoads(:,n) = targetLoad;

        % Store internal load vector
        Pest(:,n) = P;

        fprintf('\n');
        fprintf('\n      Completed.') 
        fprintf('\n      Total iterations:       %i',iter);
        fprintf('\n      Final norm(R):          %.2d',normR);
        fprintf('\n      Max. Displacement:      %.3f mm = %.2f%% of original length',max(full(Q)),max(full(Q))/Lx*100);
        fprintf('\n\n');
        
        if newAttempt ~= 0
            % Reset counter and return load increment to initial definition
            newAttempt = 0;
            DeltaP = (1/numLoadSteps).*(finalLoad - initialLoad);
        end

    end
    
    % if program failed, try an improved method
    if FAIL == 1
        if GRAD == 0 
            GRAD = 1; 
            fprintf('\n   SWITCHING TO SECANT STIFFNESS METHOD')
        
        elseif GRAD == 1 || (GRAD == 3 && newAttempt >= maxNewAttempts)
            GRAD = 2; 
            fprintf('\n   SWITCHING TO FINITE DIFFERENCE TANGENT STIFFNESS METHOD')
            n0 = n;
            
        elseif GRAD == 2 && n0 ~= n
            n0 = n;
            NEW_TANG = 0; 
            fprintf('\n   CALCULATING NEW TANGENT STIFFNESS')
            
        elseif GRAD == 2 && n0 == n
            FAIL = 2;
            fprintf('\n   FINITE DIFFERENCE TANGENT STIFFNESS METHOD INSUFFICIENT')

        elseif GRAD == 3
            % Half the load step increment for a re-try
            DeltaP = 0.5*DeltaP;
            newAttempt = newAttempt + 1;
            
        end
        n = n - 1;
    end
    
    % if program fails on tangent stiffness method, exit program
    if FAIL == 2
        fprintf('\n');
        fprintf('\n   Exiting Process Prematurely...');
        fprintf('\n      Successful Load Steps: %i',n-1);
        fprintf('\n      Percentage Complete:   %.2f%%',(n-1)/numLoadSteps*100);
        fprintf('\n');
%         outputExit;    % dump info if program failed
        break;      % exit process if last load step failed 
    else
    
    % Prepare for next load step
    n = n + 1;
    FAIL = 0;
    
    end
end

% Only consider part of solution that converged and satisfied equilibrium
if n < numLoadSteps
    numLoadSteps = n-1; % new max load
end
storeDisps = real(storeDisps(:,1:numLoadSteps));
storeLoads = storeLoads(:,1:numLoadSteps);

fprintf('\n   NON-LINEAR ELASTICITY PROBLEM COMPLETED.');
fprintf('\n');
elapsedTime1 = toc;
fprintf('\n   Time to Complete:  %.3f seconds',elapsedTime1)
fprintf('\n');

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
elseif elementType == 2
    [F,~,~,~,PK2,tau,sig,PK1,~] = postNonlin_3D(input,storeDisps,storeLoads,elementsToAnalyse);
end
elapsedTime2 = toc;
fprintf('\n   Completed!')
fprintf('\n   Time to Complete:  %.3f seconds',elapsedTime2-elapsedTime1)
fprintf('\n');

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
set(fig,'Name','Shear Stress vs Stretch (1st Piola Kirchhoff)','NumberTitle','off')
plot(storeF(:,6),storePK1(:,6),'-o');
xlabel('\gamma')
ylabel('P_{xy} [MPa]')

fig = figure;
set(fig,'Name','Lateral Stretch','NumberTitle','off')
plot(storeF(:,1),storeF(:,2),'-o');
xlabel('\lambda_{1}')
ylabel('\lambda_{2}')

fprintf('\n   Max lambda1 = %.2f',max(storeF(:,1)));
fprintf('\n   Max gamma = %.2f',max(storeF(:,6)));


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
if CASE == 3 || CASE == 91 || CASE == 92 || CASE == 93
    
    % Determine which nodes to find 
    if CASE == 91 || CASE == 92 
        endNodes = find(input.ND(:,2)==Lx);
        stress_area = Ly*Lz;
        height = Lx;
        if CASE == 91; dofI = 2; end
        if CASE == 92; dofI = 3; end

    elseif CASE == 93 || CASE == 94 || CASE == 3
        endNodes = find(input.ND(:,3)==Ly);
        stress_area = Lx*Lz;
        height = Ly;
        if CASE == 93 || CASE == 3; dofI = 1; end
        if CASE == 94; dofI = 3; end

    elseif CASE == 95 || CASE == 96
        endNodes = find(input.ND(:,4)==Lz);
        stress_area = Lx*Ly;
        height = Lz;
        if CASE == 95; dofI = 1; end
        if CASE == 96; dofI = 2; end

    end
    
    endDofX = dofPerNode.*endNodes-(dofPerNode-dofI).*ones(numel(endNodes),1);

    endStress = zeros(numLoadSteps,1); endStretch = zeros(numLoadSteps,1);
    for load = 1:numLoadSteps
        endForce = sum(storeLoads(endDofX,load));
        endStress(load) = (endForce) / (stress_area);    % tau = F/A

        endDisp = mean(storeDisps(endDofX,load));
        endStretch(load) = (endDisp) / (height);    % gam = Q/h
    end
    
    results_XX = [endStretch, endStress];

    fig = figure;
    set(fig,'Name','Calculated End Stress v Stretch');
    hold on;
    plot(endStretch,endStress);
    xlabel('Amount of Shear (gamma) = Average Displacement / Initial Height')
    ylabel('Shear Stress (PK1) = Sum of Force / Initial Area')

end



hold off

%% STRESS COUNTOURS
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


%% CONCLUDE THE SCRIPT
fprintf('\n');
fprintf('\n   ---------- PROGRAM COMPLETED ----------');
fprintf('\n');
