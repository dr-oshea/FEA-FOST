function [storeLoads, storeDisps, Pest, numLoadSteps] = nonlinear_solution_NFEA(input, elementType, dofPerNode, SolverOptions, numLoadSteps, totalDOFs, maxInc, initialLoad, finalLoad, DeltaP, K0, freeDOFs, nodesPerElem, PROG, Lx, Ly, Lz)

GRAD = SolverOptions.GRAD;       % gradient method used (3 = tangent stiffness)
BROYDEN = SolverOptions.BROYDEN;    % enable Broyden's quasi-Newton method (not relevant)
PROCMOD = SolverOptions.PROCMOD;    % enable process modification (not relevant)
QUASI_ITER = SolverOptions.QUASI_ITER; % number of iterations before new tangent stiffness calculation
TOL_equilibrium = SolverOptions.TOL_equilibrium; % tolerance on equilibrium 
TOL_procmod = SolverOptions.TOL_procmod;     % tolerance on process modifications (not relevant)
maxMod = SolverOptions.maxMod;    % maximum number of process "modifications" before exiting
maxIter = SolverOptions.maxIter;  % maximum number of iterations before exiting 

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
        if SolverOptions.MIXED == 0
            K = Assemble_Tangent_Stiffness(Q,input,elementType,nodesPerElem,dofPerNode,PROG);
        else
            [K,~,~] = Assemble_Tangent_Stiffness_MIXED(Q, input, elementType, 8, 3 ,PROG, 0, 0, 0);
        end

    end

    % Output information regarding stiffness matrix in use:
    fprintf('\n      Information regarding stiffness matrix:');
    fprintf('\n         Determinant   = %.2e',det(K(freeDOFs,freeDOFs)));
    fprintf('\n         Condition No. = %.2e',condest(K(freeDOFs,freeDOFs)));
    fprintf('\n         Symmetry      = %.2e',norm(full(K-K.')));
    fprintf('\n');
    
    Kp = K; % Store previous stiffess matrix
   
    % Perform nonlinear iterative process to reach a solution...
    [Q,P,FAIL,iter,normR] = NonlinearIteration(input,Q,P,targetLoad,K,freeDOFs,totalDOFs,nodesPerElem,TOL_equilibrium,maxIter,maxMod,QUASI_ITER,PROCMOD,FAIL,PROG,SolverOptions);  
   
    
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


end