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
elementType = 4;           % Type of element: 1 = 2D quad, 2 = 3D Brick, 3 = 2d annulus, 4 = truncated spheroid

% Geometry
Lx = 2*pi;    % [mm]    % maximum angle (should be 2*Pi)
Ly = 0.1772;       % [mm]    % thickness of annulus
Lz = 47/120;       % [mm]    % inner radius
thetaLayer1 = 0;     % [deg]
thetaLayer2 = 0;     % [deg]
numLayers = 1;        % theta1 and theta2 give orientation of material for case of multiple layers (layers > 1)

% Discretisation
meshDivisions = 2*4;   %element divisions around theta (must be a multiple of 4) %2*4   1*4
thicknessDivisions = 5; %element divisions along lambda (wall thickness) %2
heightDivisions = 3;   %element divisions along mu (height of wall) %3

eleAspectRatio = 1*4/meshDivisions; %!!CURRENTLY UNUSED!!       % change numerator for element divisions along radius [keeping aspect ratio = 1 for now] eleAspectRatio = 4/meshDivisions

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
    use_dispersion = 0; % 1 = on, 0 = off

    if constitutiveModel == 1   % Myocardium model (2022)
        L00 = 10;  
        LFF = 717;
        LSS = 215;
        LNN = 95;
        M00 = 3;
        MFS = 0;
        if use_dispersion == 1
            kappaF = 0.1738; %0.1738
            kappaS = 0.1179; %0.1738
        else
            kappaF = 0.0; %0.000
            kappaS = 0.0; %0.000
        end
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
        if use_dispersion == 1
            kappaF = 0.08; 
            kappaS = 0.09;
        else
            kappaF = 0;
            kappaS = 0;
        end
        input.model = [constitutiveModel, a, b, af, bf, as, bs, afs, bfs, kappaF, kappaS];

    end
end
PLANAR = 2; % 1 = Plane stress, 2 = Plane strain;

% Specified final loading (Force OR Stress used, depending on problem)
numLoadSteps = 10;      % Load steps (increasing takes longer, though more stable solution process) - number of divisions (actual load steps = N+1) %2
finalForce = 5.0;      % [Newtons]
finalStress = 15;    % [should match units of material parameters, kPa, MPa, etc..] %0.01

% Solution procedure options
GRAD = 3;       % gradient method used (3 = tangent stiffness)
BROYDEN = 0;    % enable Broyden's quasi-Newton method (not relevant)
PROCMOD = 1;    % enable process modification (not relevant)
QUASI_ITER = 1; % number of iterations before new tangent stiffness calculation
TOL_equilibrium = 10^-3; % tolerance on equilibrium 
TOL_procmod = 10^-2;     % tolerance on process modifications (not relevant)
maxMod = 10;    % maximum number of process "modifications" before exiting
maxIter = 11;  % maximum number of iterations before exiting %was originally 1000
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
input.ND = input_nodal(meshDivisions,eleAspectRatio,Lx,Ly,Lz,elementType,thicknessDivisions,heightDivisions);

% Array storing all element connectivity
[input.EL, input.FIBRES] = input_element(input.ND,meshDivisions,eleAspectRatio,numLayers,Lx,Ly,Lz,thetaLayer1,thetaLayer2,elementType,thicknessDivisions,heightDivisions,use_dispersion);

% Array storing all material properties and material tensors
% MATERIAL PROPERTIES ARE HARD CODED IN in FILE Get_S_CC.m
% [input.C0,input.L0,input.M0] = input_material(MatProps,PLANAR,type);

% Array storing all boundary conditions
input.CON  = input_constraint(input.ND,Lx,Ly,Lz,elementType,CCASE,meshDivisions,eleAspectRatio,thicknessDivisions,heightDivisions);
input.LOAD = input_load(input.ND,input.EL,Lx,Ly,Lz,meshDivisions,elementType,finalForce,finalStress,FCASE,meshDivisions,eleAspectRatio,thicknessDivisions,heightDivisions);

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

%fig = figure;
% set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
plotDiscretisationShaded(input) %plots blank discretisation with shading so its even nicer to look at
% 
% fig = figure;
% set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
% plotDiscretisationBlank(input,0,elementType) %just plots the 3d mesh without any shading
% 
fig = figure;
set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
plotDiscretisation(input,1,elementType)
% 
% fig = figure;
% set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
% plotConstraints(input,0,elementType)
% 
% fig = figure;
% set(fig,'Name','FEM Fibre Dispersion','NumberTitle','off');
% plotFibre(input,0,elementType) %only for element type 4
% 
% fig = figure;
% set(fig,'Name','FEM Sheet Orientation','NumberTitle','off');
% plotSheet(input,0,elementType) %only for element type 4
% 
% fig = figure;
% set(fig,'Name','FEM Sheet Normal Orientation','NumberTitle','off');
% plotSheetNormal(input,0,elementType) %only for element type 4
% 
% fig = figure;
% set(fig,'Name','FEM Loading','NumberTitle','off');
% plotLoad(input,0,elementType) %only for element type 4
% 
% fig = figure;
% set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
% plotLoadShaded(input)
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
maxNewAttempts = 10; %i think this the number of times loadstep is haved before switching methods

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
            fprintf('\n   HALVING THE LOAD STEP AND TRYING AGAIN')
            
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

%% Save everything you need
s = 'saved_location.txt';
dPath = which (s);
cd (dPath(1:end-length(s)));

save("input.mat", "input")
save("displacements.mat", "storeDisps")
save("loads.mat", "storeLoads")

clear


% Call postprocessing function.
post_proc = 1;
if post_proc == 1; post_processing_kiru; end



%% Plot stretch vs stress, calculated using force v displacement


%% CONCLUDE THE SCRIPT
fprintf('\n');
fprintf('\n   ---------- PROGRAM COMPLETED ----------');
fprintf('\n');
