%% FEA_nonlin2D

% This script runs the nonlinear finite element program where the strain
% energy function and corresponding constitutive law is described by the 
% work submitted by A/Prof. Mario Attard and Daniel O'Shea in 2018. The 
% program uses the Initial or Secant Gradient Method to iteratively solve a
% nonlinear hyperelasticity problem in 2D/3D subject to given applied nodal 
% loads and constraints

% Author:   Daniel O'Shea   < d.oshea@unsw.edu.au >
% Created:  16 March 2018

%%                             Preliminary

% Clear the workspace
clc
clear
close all
global nset cset

% Select Case to Study:
% 1 = Pinned bottom edge, tension top edge Y- direction
% 2 = Uniaxial tension in X-Direction
% 3 = Ciarletta shear test
% 4 = Simple shear test
% 5 = Equibiaxial tension
% 6 = Uniaxial Tension 1 Element
% 7 = Pure Shear (Biaxial Testing)
CASE = 3;
MAT = 3;    % 1 = specify isotropic, 2 = specify orthotropic, 3 = known sets


fprintf('\n   Commencing Nonlinear FE Analysis Program');
fprintf('\n');

% Start Timer:
tic

%%                           Define the problem

% This section defines the geometry of the continuum, the discretisation,
% the material properties, and the boundary conditions

% Dimension of problem: Type of Element used
type = 1;           % Type of element: 1 = 2D quad, 2 = 3D Brick

% Geometry
xdim = 10;       % [mm]
ydim = 5;       % [mm]
zdim = 1;       % [mm]
theta1 = 0;     % [deg]
theta2 = 0;     % [deg]
layers = 1;

% Discretisation
mesh_div = 2;       % element divisions along x-dimension (2D), z-dim (3D)
aspect = 1;         % aspect ratio of elements

% Isotropic Material Properties
if MAT == 1
    E = 1.434;
    v = 0.495;
    G = E/(2*(1+v));
    MatProps = [E E E v v v G G G];
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
    MatProps = [E11 E22 E33 v23 v13 v12 G23 G13 G12];
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
    MatProps = Ciarletta2009;
    
end
PLANAR = 2; % 1 = Plane stress, 2 = Plane strain;

% Specified final loading
N = 500;        % Load steps
Force = 0.0;      % [Newtons]
Stress = 0.0012;    % [Newtons / mm / mm]

% Finite strain indices
n1 = 0.750;
n2 = -0.703;
n3 = 2.628;
% nset = [n1; n2; n3];
nset = [4.3629];

% Mu additional constants
c2 = 0.0;
% cset = [1-c2-c3; c2; c3];
cset = [1];

% Solution options
GRAD = 1;       % gradient method used
BROYDEN = 0;    % enable Broyden's quasi-Newton method
PROCMOD = 1;    % enable process modification
alpha1 = 10^-3; % tolerance on equilibrium
alpha2 = 10^-2; % tolerance on process modifications
maxMOD = 10;    % maximum number of modifications
maxITER = 300; % maximum number of iterations
PLOT = 0;
PROG = 0;       % display progress bars

%%                           Generate input arrays
fprintf('\n   Determining Input Information...');

% Boundary conditions
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
    CCASE = 5;
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
end
fprintf('\n');

% Array storing all nodal values
input.ND = input_nodal(mesh_div,aspect,xdim,ydim,zdim,type);

% Array storing all element connectivity
input.EL = input_element(input.ND,mesh_div,layers,xdim,ydim,zdim,theta1,theta2,type);

% Array storing all material properties and material tensors
[input.C0,input.L0,input.M0] = input_material(MatProps,PLANAR,type);

% Array storing all boundary conditions
input.CON  = input_constraint(input.ND,xdim,ydim,zdim,type,CCASE);
input.LOAD = input_load(input.ND,input.EL,xdim,ydim,zdim,mesh_div,type,Force,Stress,FCASE);

% Array storing element information
[dofN, ndsE, strC] = input_info(type);

% Additional info
numN = size(input.ND,1);
dofT = numN*dofN;
numE = size(input.EL,1);


%% Display details of defined elasticity problem
fprintf('\n   Problem Generated.');
fprintf('\n');
fprintf('\n   Summary of Problem:');
fprintf('\n      Body width  = %.2f mm',xdim);
fprintf('\n      Body height = %.2f mm',ydim);
fprintf('\n      Constriants:');
if CCASE == 1
    fprintf('\n         Bottom edge pinned.');
elseif CCASE == 2
    fprintf('\n         Centre node pinned.');
end
fprintf('\n      Loading:');
if FCASE == 1
    fprintf('\n         Top edge uniform load = %.2f N/mm2',Stress);
elseif FCASE == 2
    fprintf('\n         Top edge uniform load    =  %.2f N/mm2',Stress);
    fprintf('\n         Bottom edge uniform load = -%.2f N/mm2',Stress);
end
fprintf('\n');
fprintf('\n   Summary of Discretisation:');
if type == 1
    fprintf('\n      2D Bilinear Quadrilateral Elements');
elseif type == 2
    fprintf('\n      3D Hexahedral Elements');
end
fprintf('\n      Number of Nodes:          %i',size(input.ND,1));
fprintf('\n      Number of Elements:       %i',size(input.EL,1));
fprintf('\n      Total Degrees of Freedom: %i',dofT);
fprintf('\n');
fprintf('\n   Input Material Properties:');
fprintf('\n      E11 = %.4f MPa',MatProps(1));
fprintf('\n      E22 = %.4f MPa',MatProps(2));
fprintf('\n      E33 = %.4f MPa',MatProps(3));
fprintf('\n      v23 = %.4f    ',MatProps(4));
fprintf('\n      v13 = %.4f    ',MatProps(5));
fprintf('\n      v12 = %.4f    ',MatProps(6));
fprintf('\n      G23 = %.4f MPa',MatProps(7));
fprintf('\n      G13 = %.4f MPa',MatProps(8));
fprintf('\n      G12 = %.4f MPa',MatProps(9));
fprintf('\n');


%%                           Plot of Discretisation

% This section produces a plot of the discretised mesh and applied boundary
% condtitions

fig = figure;
set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
plotDiscretisation(input,0,type)

%%                         Determine Initial Gradient

% This section determines the initial gradient of the elasticity problem.
% This is equivalent to the stiffness matrix of each element in the
% linear infinitesimal strain region

fprintf('\n   Determining Initial Stiffness Matrix...');

K0 = Assemble_Initial_Stiffness(input,type,ndsE,dofN);

fprintf('\n   Completed.');
fprintf('\n');
fprintf('\n   Summary of Initial Stiffness Matrix:');
fprintf('\n      Total size:             %i x %i',dofN*size(input.ND,1),dofN*size(input.ND,1));
fprintf('\n      Total entries:          %i',(dofN*size(input.ND,1))^2);
fprintf('\n      Total non-zero entries: %i',nnz(K0));
fprintf('\n');


%%                   Define Incremental Load Vector

% This section generates the final load vector based on input conditions,
% before using this to determine the global incremental load vector

fprintf('\n   Constructing incremental load vector...');

% Number of load steps
fprintf('\n      Total load steps:     %i',N);

% Final Load Vector
PN = Final_Load_Vector(input);
fprintf('\n      Final load vector constructed');
fprintf('\n      Maximum Load:         %.3f N',max(max(input.LOAD(:,2:end))));

% Initial Load Vector
P0 = 0.*PN;

% Incremental load vector
DeltaP = (1/N).*(PN - P0);
fprintf('\n      Incremental load vector determined');
fprintf('\n');
maxInc = max(max(input.LOAD(:,2:end)))/N;

%%                     Commence Iterative Process

% Determine vector of unconstrained degrees of freedom
L = Apply_Constraints(input);


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
    fprintf('\n   "Tangent Stiffness Method" \n');    
end

% Allocate initial displacement and load vectors (zero vectors)
Qsol = zeros(dofT,N);
Q0   = zeros(dofT,1);
Psol = zeros(dofT,N);
Pest = zeros(dofT,N);

% Solution tolerance
fprintf('\n   Solution Tolerance: %.0d',alpha1)
fprintf('\n')

% Iterative Process:
maxTar = maxInc;

fprintf('\n   Process Modification:');
fprintf('\n      Eigenvalue Tolerance:  %.0d',alpha2);
fprintf('\n      Maximum modifications: %i',maxMOD);
fprintf('\n');

Beta = zeros(maxMOD,1);
Phi = zeros(dofT,maxMOD);
NEW_TANG = 0;
FAIL = 0;
Kp = zeros(dofT,dofT);

% For each load step
n = 1;
while n <= N
    fprintf('\n\n   COMMENCING LOAD STEP #%i',n);
    
    % Initial Guess = previous known solution
    if n == 1
        Q = Q0;
        P = P0;
    else
        Q = Qsol(:,n-1);     % previous known vector
        P = Psol(:,n-1)';    % previous known vector
    end
    
    % Targeted Load Vector = previous known + deltaP
    Pn = P + DeltaP;
    
    % Initial Stiffness
    if GRAD == 0 || (GRAD == 2 && n == 1)
        K = K0;
    
    % Secant Stiffness
    elseif GRAD == 1
        fprintf('\n      Assembling secant stiffness matrix...')
        K = Assemble_Secant_Stiffness(Q,input,type,ndsE,dofN);
    
    % Finite Difference Tangent Stiffness
    elseif GRAD == 2 && NEW_TANG == 0
        fprintf('\n      Assembling tangent stiffness matrix...')
        hwb    = waitbar(0,'Assembly process ...');
        set(hwb,'Name','Approximating Tangent Stiffness Matrix');
        indj = 1;
        eta = 10^-3;
        for j = 1:dofT
            if j>(indj*10)
                waitbar(j/dofT,hwb,[num2str(floor(100*j/dofT)) ' % assembly ...']);
                indj=indj+1; 
            end
            if j==dofT
                waitbar(j/dofT,hwb,' assembly finished ');
                disp(['     K assembly: ' num2str(i) ' dofs - complete']);  
            end
            signQ = 1;
            if Q(j) < 0; signQ = -1; end
            h = sqrt(eta)*max(abs(Q(j)),1)*signQ;
            Qj = Q(j);
            Q(j) = Q(j) + h;
            h = Q(j) - Qj;
            Pj = SolveP(Q,input,ndsE);
            for i = 1:dofT
                K(i,j) = (Pj(i)-P(i))/h;
            end
            Q(j) = Qj;
        end
        NEW_TANG = 1;
        close(hwb);
        
    end
    fprintf('\n      Information regarding stiffness matrix:');
    fprintf('\n         Determinant   = %.2e',det(K(L,L)));
%     fprintf('\n         Max. Diagonal = %.2d',max(diag(full(K))));
%     fprintf('\n         Min. Diagonal = %.2d',min(diag(full(K))));
    fprintf('\n         Condition No. = %.2e',condest(K(L,L)));
    fprintf('\n         Symmetry      = %.2e',norm(full(K-K.')));
    fprintf('\n');
    Kp = K; % Store previous stiffess matrix
    
    [Q,P,FAIL,iter,normR] = NonlinearIteration(input,Q,P,Pn,K,L,dofT,ndsE,alpha1,maxITER,maxMOD,PROCMOD,FAIL);  
    
    % if program completed current load step
    if FAIL == 0
        % Store converged displacement solution
        Qsol(:,n) = Q;

        % Store target load vector
        Psol(:,n) = Pn;

        % Store calculated load vector
        Pest(:,n) = P;

        fprintf('\n');
        fprintf('\n      Completed.') 
        fprintf('\n      Total iterations:       %i',iter);
        fprintf('\n      Final norm(R):          %.2d',normR);
        fprintf('\n      Max. Displacement:      %.3f mm = %.2f%% of original length',max(full(Q)),max(full(Q))/xdim*100);
        fprintf('\n\n');
    end
    
    % if program failed, try an improved method
    if FAIL == 1
        if GRAD == 0 
            GRAD = 1; 
            fprintf('\n   SWITCHING TO SECANT STIFFNESS METHOD')
        
        elseif GRAD == 1 
            GRAD = 2; 
            fprintf('\n   SWITCHING TO TANGENT STIFFNESS METHOD')
            n0 = n;
            
        elseif GRAD == 2 && n0 ~= n
            n0 = n;
            NEW_TANG = 0; 
            fprintf('\n   CALCULATING NEW TANGENT STIFFNESS')
            
        elseif GRAD == 2 && n0 == n
            FAIL = 2;
            fprintf('\n   TANGENT STIFFNESS METHOD INSUFFICIENT')
            
        end
        n = n - 1;
    end
    
    % if program fails on tangent stiffness method, exit program
    if FAIL == 2
        fprintf('\n');
        fprintf('\n   Exiting Process Prematurely...');
        fprintf('\n      Successful Load Steps: %i',n-1);
        fprintf('\n      Percentage Complete:   %.2f%%',(n-1)/N*100);
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
if n < N
    N = n-1; % new max load
end
Qsol = real(Qsol(:,1:N));
Psol = Psol(:,1:N);

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
ele_analyse = input.EL(:,1);    % all elements

% Run post processing (returns deformation and stress vectors)
if type == 1
    [F,~,~,~,~,tau,P,~] = postNonlin_2D(input,Qsol,Psol,ele_analyse);
elseif type == 2
    [F,~,~,~,~,tau,P,~] = postNonlin_3D(input,Qsol,Psol,ele_analyse);
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
def = [zeros(size(input.ND,1),1) reshape(Qsol(:,end),dofN,size(input.ND,1))'] + coords0;
xmin = min([min(input.ND(:,2)),min(def(:,2))]);
xmax = max([max(input.ND(:,2)),max(def(:,2))]);
ymin = min([min(input.ND(:,3)),min(def(:,3))]);
ymax = max([max(input.ND(:,3)),max(def(:,3))]);
if type == 2
    zmin = min([min(input.ND(:,4)),min(def(:,4))]);
    zmax = max([max(input.ND(:,4)),max(def(:,4))]);
end

% Plot Deformed Shapes
for i = N:-max(factor(N)):1
    
    Q = Qsol(:,i);
    
    Qmap = reshape(Q,dofN,size(input.ND,1))';
    
    input.ND(:,2:end) = coords0(:,2:end) + Qmap;
    
    fig = figure;
    str = sprintf('FEM Solution: Load Step #%i',i);
    set(fig,'Name',str,'NumberTitle','off');
    plotDiscretisation(input,0,type);
    xlim([xmin xmax])
    ylim([ymin ymax])
    if type == 2
        zlim([zmin zmax]);
    end
    
end
input.ND = coords0; % reset the input data to original coordinates

%% plot Stress v Stretch

fprintf('\n   Plotting stress vs stretch curve...');
fprintf('\n');

numGrad = dofN * dofN;
xmin = min(input.ND(:,2));
ymin = min(input.ND(:,3));
spacing = (xdim-xmin)/mesh_div;

% node to plot stress v stretch at
if type == 1
    node = find(input.ND(:,2)==xdim/2 & input.ND(:,3)==ydim/2);
elseif type == 2
    node = find(input.ND(:,2)==xdim/2 & input.ND(:,3)==ydim/2 & input.ND(:,4)==zdim/2);
end

lamx = zeros(N,1); tauxx = zeros(N,1); Pxx = zeros(N,1);
lamy = zeros(N,1); tauyy = zeros(N,1); Pxy = zeros(N,1);
gamxy = zeros(N,1);
for i = 1 : N
    
    % Get elements connected to specified node
    el1 = find(input.EL(:,1+2) == node);
    el2 = find(input.EL(:,2+2) == node);
    el3 = find(input.EL(:,3+2) == node);
    el4 = find(input.EL(:,4+2) == node);
    
    % Get F11 from each element (lambda1)
    v1 = F(el1,numGrad*(1-1)+1,i);
    v2 = F(el2,numGrad*(2-1)+1,i);
    v3 = F(el3,numGrad*(3-1)+1,i);
    v4 = F(el4,numGrad*(4-1)+1,i);
    lamx(i) = mean([v1 v2 v3 v4]);

    % Get F22 from each element (lambda2)
    v1 = F(el1,numGrad*(1-1)+2,i);
    v2 = F(el2,numGrad*(2-1)+2,i);
    v3 = F(el3,numGrad*(3-1)+2,i);
    v4 = F(el4,numGrad*(4-1)+2,i);
    lamy(i) = mean([v1 v2 v3 v4]);
    
    % Get F12 from each element (gamma for shear)
    v1 = F(el1,numGrad*(1-1)+3,i);
    v2 = F(el2,numGrad*(2-1)+3,i);
    v3 = F(el3,numGrad*(3-1)+3,i);
    v4 = F(el4,numGrad*(4-1)+3,i);
    gamxy(i) = mean([v1 v2 v3 v4]);
    
    % Get tau11 from each element
    v1 = tau(el1,numGrad*(1-1)+1,i);
    v2 = tau(el2,numGrad*(2-1)+1,i);
    v3 = tau(el3,numGrad*(3-1)+1,i);
    v4 = tau(el4,numGrad*(4-1)+1,i);
    tauxx(i) = mean([v1 v2 v3 v4]);
    
    % Get tau22 from each element
    v1 = tau(el1,numGrad*(1-1)+2,i);
    v2 = tau(el2,numGrad*(2-1)+2,i);
    v3 = tau(el3,numGrad*(3-1)+2,i);
    v4 = tau(el4,numGrad*(4-1)+2,i);
    tauyy(i) = mean([v1 v2 v3 v4]);
    
    % Get P11 from each element
    v1 = P(el1,numGrad*(1-1)+1,i);
    v2 = P(el2,numGrad*(2-1)+1,i);
    v3 = P(el3,numGrad*(3-1)+1,i);
    v4 = P(el4,numGrad*(4-1)+1,i);
    Pxx(i) = mean([v1 v2 v3 v4]);
    
    % Get P12 from each element
    v1 = P(el1,numGrad*(1-1)+3,i);
    v2 = P(el2,numGrad*(2-1)+3,i);
    v3 = P(el3,numGrad*(3-1)+3,i);
    v4 = P(el4,numGrad*(4-1)+3,i);
    Pxy(i) = mean([v1 v2 v3 v4]);
           
end    
    
fig = figure;
set(fig,'Name','Stress vs Stretch (XX Component)','NumberTitle','off')
plot(lamx,tauxx,'-o');
xlabel('\lambda_{1}')
ylabel('\tau_{xx}  [MPa]')

fig = figure;
set(fig,'Name','Stress vs Stretch (YY Component)','NumberTitle','off')
plot(lamx,tauyy,'-o');
xlabel('\lambda_{1}')
ylabel('\tau_{yy}  [MPa]')

fig = figure;
set(fig,'Name','Tensile Stress vs Stretch (1st Piola Kirchhoff)','NumberTitle','off')
plot(lamx,Pxx,'-o');
xlabel('\lambda_{1}')
ylabel('P_{xx} [MPa]')

fig = figure;
set(fig,'Name','Shear Stress vs Stretch (1st Piola Kirchhoff)','NumberTitle','off')
plot(gamxy,Pxy,'-o');
xlabel('\gamma')
ylabel('P_{xy} [MPa]')

fig = figure;
set(fig,'Name','Lateral Stretch','NumberTitle','off')
plot(lamx,lamy,'-o');
xlabel('\lambda_{1}')
ylabel('\lambda_{2}')

fprintf('\n   Max lambda1 = %.2f',max(lamx));
fprintf('\n   Max gamma = %.2f',max(gamxy));

%% Plot stretch vs stress, calculated using force v displacement

% Uniaxial Tension Plotting
if CASE == 2 || CASE == 6 || CASE == 5
    endnodes = find(input.ND(:,2)==xdim);
    enddofx = dofN.*endnodes-(dofN-1).*ones(numel(endnodes),1);

    EndStress = zeros(N,1); EndStretch = zeros(N,1);
    for i = 1:N
        EndForce = sum(Psol(enddofx,i));
        EndStress(i) = (2*EndForce) / (2*ydim*zdim);    % tau = F/A

        EndDisp = mean(Qsol(enddofx,i));
        EndStretch(i) = (2*EndDisp) / (xdim) + 1;    % lam = Q/L + 1
    end

    fig = figure;
    set(fig,'Name','Calculated End Stress v Stretch');
    plot(EndStretch,EndStress);
    xlabel('Stretch = 1 + 2*ave(Qx)/xdim')
    ylabel('P11 = 2*sum(Fx)/2*ydim')
end

% Ciarletta-type Shear Experiment
if CASE == 3
    endnodes = find(input.ND(:,3)==ydim);
    enddofx = dofN.*endnodes-(dofN-1).*ones(numel(endnodes),1);

    EndStress = zeros(N,1); EndStretch = zeros(N,1);
    for i = 1:N
        EndForce = sum(Psol(enddofx,i));
        EndStress(i) = (EndForce) / (xdim*zdim);    % tau = F/A

        EndDisp = mean(Qsol(enddofx,i));
        EndStretch(i) = (EndDisp) / (ydim);    % gam = Q/h
    end

    fig = figure;
    set(fig,'Name','Calculated End Stress v Stretch');
    hold on;
    plot(EndStretch,EndStress);
    xlabel('gamma = ave(Qx)/ydim')
    ylabel('P12 = sum(Fx)/(xdim*zdim)')
    xlim([0 0.2])
    ylim([0 0.045])
end

hold off

%% STRESS COUNTOURS
if type == 1 && PLOT ~= 0
Xloc = reshape(input.ND(:,2),[mesh_div+1,mesh_div+1]);
Yloc = reshape(input.ND(:,3),[mesh_div+1,mesh_div+1]);

% CALCULATE ALL ELEMENT VALUES OF COMPONENTS
ele_analyse = input.EL(:,1);
[F,C,E0,En,S,tau,P,q] = postNonlin_2D(input,Qsol,Psol,ele_analyse);

% Specify which load steps to plot contours at
load_analyse = [1; 2; floor(N/2); (N-10:1:N)'];
scrsz = get(groot,'ScreenSize');

for i = 1:numel(load_analyse)
    ld = load_analyse(i);
    lam1 = zeros(numN,1);
    lam2 = zeros(numN,1);
    tau11 = zeros(numN,1);
    tau22 = zeros(numN,1);
    tau12 = zeros(numN,1);
    E011 = zeros(numN,1);
    E022 = zeros(numN,1);
    En11 = zeros(numN,1);
    En22 = zeros(numN,1);
    
    for a = 1:numN
        node = input.ND(a,1);
    
        % Get elements connected to specified node
        el1 = find(input.EL(:,1+2) == node);
        el2 = find(input.EL(:,2+2) == node);
        el3 = find(input.EL(:,3+2) == node);
        el4 = find(input.EL(:,4+2) == node);
    
        % Get F11 from each element (lambda1)
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = F(el1,numGrad*(1-1)+1,ld); end
        if isempty(el2)==0; v2 = F(el2,numGrad*(2-1)+1,ld); end
        if isempty(el3)==0; v3 = F(el3,numGrad*(3-1)+1,ld); end
        if isempty(el4)==0; v4 = F(el4,numGrad*(4-1)+1,ld); end
        lam1(a) = mean([v1 v2 v3 v4]);
        
        % Get F22 from each element (lambda2)
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = F(el1,numGrad*(1-1)+2,ld); end
        if isempty(el2)==0; v2 = F(el2,numGrad*(2-1)+2,ld); end
        if isempty(el3)==0; v3 = F(el3,numGrad*(3-1)+2,ld); end
        if isempty(el4)==0; v4 = F(el4,numGrad*(4-1)+2,ld); end
        lam2(a) = mean([v1 v2 v3 v4]);

        % Get E011 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = E0(el1,numGrad*(1-1)+1,ld); end
        if isempty(el2)==0; v2 = E0(el2,numGrad*(2-1)+1,ld); end
        if isempty(el3)==0; v3 = E0(el3,numGrad*(3-1)+1,ld); end
        if isempty(el4)==0; v4 = E0(el4,numGrad*(4-1)+1,ld); end
        E011(a) = mean([v1 v2 v3 v4]);
        
        % Get E022 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = E0(el1,numGrad*(1-1)+2,ld); end
        if isempty(el2)==0; v2 = E0(el2,numGrad*(2-1)+2,ld); end
        if isempty(el3)==0; v3 = E0(el3,numGrad*(3-1)+2,ld); end
        if isempty(el4)==0; v4 = E0(el4,numGrad*(4-1)+2,ld); end
        E022(a) = mean([v1 v2 v3 v4]);
        
        % Get En11 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = En(el1,numGrad*(1-1)+1,ld); end
        if isempty(el2)==0; v2 = En(el2,numGrad*(2-1)+1,ld); end
        if isempty(el3)==0; v3 = En(el3,numGrad*(3-1)+1,ld); end
        if isempty(el4)==0; v4 = En(el4,numGrad*(4-1)+1,ld); end
        En11(a) = mean([v1 v2 v3 v4]);
        
        % Get En22 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = En(el1,numGrad*(1-1)+2,ld); end
        if isempty(el2)==0; v2 = En(el2,numGrad*(2-1)+2,ld); end
        if isempty(el3)==0; v3 = En(el3,numGrad*(3-1)+2,ld); end
        if isempty(el4)==0; v4 = En(el4,numGrad*(4-1)+2,ld); end
        En22(a) = mean([v1 v2 v3 v4]);
        
        % Get tau11 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = tau(el1,numGrad*(1-1)+1,ld); end
        if isempty(el2)==0; v2 = tau(el2,numGrad*(2-1)+1,ld); end
        if isempty(el3)==0; v3 = tau(el3,numGrad*(3-1)+1,ld); end
        if isempty(el4)==0; v4 = tau(el4,numGrad*(4-1)+1,ld); end
        tau11(a) = mean([v1 v2 v3 v4]);
        
        % Get tau22 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = tau(el1,numGrad*(1-1)+2,ld); end
        if isempty(el2)==0; v2 = tau(el2,numGrad*(2-1)+2,ld); end
        if isempty(el3)==0; v3 = tau(el3,numGrad*(3-1)+2,ld); end
        if isempty(el4)==0; v4 = tau(el4,numGrad*(4-1)+2,ld); end
        tau22(a) = mean([v1 v2 v3 v4]);
        
        % Get tau12 from each element
        v1=[]; v2=[]; v3=[]; v4=[];
        if isempty(el1)==0; v1 = tau(el1,numGrad*(1-1)+3,ld); end
        if isempty(el2)==0; v2 = tau(el2,numGrad*(2-1)+3,ld); end
        if isempty(el3)==0; v3 = tau(el3,numGrad*(3-1)+3,ld); end
        if isempty(el4)==0; v4 = tau(el4,numGrad*(4-1)+3,ld); end
        tau12(a) = mean([v1 v2 v3 v4]);
        
    end
    
    fig = figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);
    str = sprintf('Contour Plots: Load Step #%i',ld);
    set(fig,'Name',str,'NumberTitle','off');

    subplot(2,4,1)
    Z = reshape(lam1,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('Lambda 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(lam1)));
    xlabel(str)
    colorbar
    
    subplot(2,4,2)
    Z = reshape(lam2,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('Lambda 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(lam2)));
    xlabel(str)
    colorbar
    
    subplot(2,4,3)
    Z = reshape(tau11,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('tau 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(tau11)));
    xlabel(str)
    colorbar
    
    subplot(2,4,4)
    Z = reshape(tau22,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('tau 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(tau22)));
    xlabel(str)
    colorbar
    
    subplot(2,4,5)
    Z = reshape(E011,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('E0 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(E011)));
    xlabel(str)
    colorbar
    
    subplot(2,4,6)
    Z = reshape(E022,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('E0 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(E022)));
    xlabel(str)
    colorbar
    
    subplot(2,4,7)
    Z = reshape(En11,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('En 11')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(En11)));
    xlabel(str)
    colorbar
    
    subplot(2,4,8)
    Z = reshape(En22,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('En 22')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(En22)));
    xlabel(str)
    colorbar
    
    figure
    Z = reshape(tau12,[mesh_div+1,mesh_div+1]);
    contourf(Xloc,Yloc,Z,'ShowText','on');
    title('tau 12')
    str = sprintf('Average Absolute Value = %.2d',mean(abs(tau12)));
    xlabel(str)
    colorbar
    
    
end
end
%% COMPARE WITH HYPERELASTIC MODEL



%% CONCLUDE THE SCRIPT
fprintf('\n');
fprintf('\n   ---------- PROGRAM COMPLETED ----------');
fprintf('\n');
