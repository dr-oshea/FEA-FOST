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
CASE = 6;   % choose boundary conditions from previously set-up cases.
MAT = 4;    % 1 = specify isotropic, 2 = specify orthotropic, 3 = known sets of properties

fprintf('\n   Commencing Nonlinear FE Analysis Program');
fprintf('\n');

% Start Timer:
tic

%%                           Define the problem

% This section defines the geometry of the continuum, the discretisation,
% the material properties, and the boundary conditions (loads and
% constraints)

% Dimension of problem: Type of Element used
elementType = 2;           % Type of element: 1 = 2D quad, 2 = 3D Brick

% Geometry
Lx = 1;       % [mm]
Ly = 1;       % [mm]
Lz = 1;       % [mm]
thetaLayer1 = 0;     % [deg]
thetaLayer2 = 0;     % [deg]
numLayers = 1;        % theta1 and theta2 give orientation of material for case of multiple layers (layers > 1)

% Discretisation
meshDivisions = 1;       % element divisions along x-dimension (2D), z-dim (3D)
eleAspectRatio = 1;     % aspect ratio of elements (1 = square/cubic)

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
        L00 = 10;  % this needs to be orders of magnitude larger than other values (enforces incompressibility)
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
finalStress = 0.5;      % [should match units of material parameters, kPa, MPa, etc..]

% Solution procedure options
SolverOptions.GRAD = 3;       % gradient method used (3 = tangent stiffness)
SolverOptions.BROYDEN = 0;    % enable Broyden's quasi-Newton method (not relevant)
SolverOptions.PROCMOD = 1;    % enable process modification (not relevant)
SolverOptions.QUASI_ITER = 5; % number of iterations before new tangent stiffness calculation
SolverOptions.TOL_equilibrium = 10^-3; % tolerance on equilibrium 
SolverOptions.TOL_procmod = 10^-2;     % tolerance on process modifications (not relevant)
SolverOptions.maxMod = 10;    % maximum number of process "modifications" before exiting
SolverOptions.maxIter = 100;  % maximum number of iterations before exiting 
SolverOptions.displayProgress = "off";

PLOT = 0;       % choose whether to suppress plotting (not relevant)
PROG = 0;       % display progress bars (this works! - set 0 to not show)

SolverOptions.lineSearchOptions.use = 1;
SolverOptions.lineSearchOptions.method = "Quadratic Expansion";
SolverOptions.lineSearchOptions.tolerance = 0.8;
SolverOptions.lineSearchOptions.maxSearch = 15;

SolverOptions.MIXED = 1;    % Use mixed formulation for incompressibility

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

% Set up "input" object: contains infor on Nodes, Elements, Loads, and
% Constraints

% Array storing all nodal values
input.ND = input_nodal(meshDivisions,eleAspectRatio,Lx,Ly,Lz,elementType);

% Array storing all element connectivity
[input.EL, input.FIBRES] = input_element(input.ND,meshDivisions,eleAspectRatio,numLayers,Lx,Ly,Lz,thetaLayer1,thetaLayer2,elementType);

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

%%      Take input, and Execute Pre-processing

[K0, DeltaP, initialLoad, finalLoad, maxInc] = pre_processing_NFEA(input, elementType, Lx, Ly, Lz, CCASE, FCASE, totalDOFs, nodesPerElem, dofPerNode, numLoadSteps, finalStress);


%%                     Commence Iterative Process

% Determine vector of unconstrained degrees of freedom
freeDOFs = Apply_Constraints(input);


%%                         Solve elasticity problem

% This section iteratively estimates the displacement solution at each
% defined load step until specified convergence limits are met

[storeLoads, storeDisps, Pest, numLoadSteps] = nonlinear_solution_NFEA(input, elementType, dofPerNode, SolverOptions,numLoadSteps, totalDOFs, maxInc, initialLoad, finalLoad, DeltaP, K0, freeDOFs, nodesPerElem, PROG, Lx, Ly, Lz);

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

%% Build functionality of saving results...


%%                      Plotting of Results

% End game - would be nice to plot changing displacement/stress contour
% plots as load step progresses (use app designer)

post_processing_NFEA(input, CASE, storeDisps, storeLoads, numLoadSteps, F, PK2, tau, sig, PK1, dofPerNode, nodesPerElem, elementType, Lx, Ly, Lz, meshDivisions);


%% CONCLUDE THE SCRIPT
fprintf('\n');
fprintf('\n   ---------- PROGRAM COMPLETED ----------');
fprintf('\n');
