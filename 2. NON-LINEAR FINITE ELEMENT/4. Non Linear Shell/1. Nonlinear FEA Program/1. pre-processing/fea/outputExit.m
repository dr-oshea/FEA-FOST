
%% DUMP OUT SOME RESULTS FOR A GIVEN LOAD STEP

% Output stiffness matrix, the diagonal values, max and min, the
% determinant, number of terms, number of nonzero terms

% Maybe begin iterations again starting from this load step but outputting
% all the load vectors, the C tensors etc....

fprintf('\n    COMMENCING DETAILED ANALYSIS OF FAILED LOAD STEP');

% run iterative process for current load step again, outputting the C
% tensor, E tensor, for a given element at a given point (centre of an
% element)???

fprintf('\n    Restarting the iterative process, outputting all relevant variables');

% Initial Guess = previous known solution
if n == 1
    Q = Q0;
    P = P0;
else
    Q = Qsol(:,n-1);     % previous known vector
    P = Psol(:,n-1)';    % previous known vector
end
fprintf('\n  Previous Dispalcement Solution:\n');
disp(Q)

% Targeted Load Vector = previous known + deltaP
Pn = P + DeltaP;
fprintf('\n  Targeted Load Vector:\n');
disp(Pn')

% Initial Residual = Target Load - Previous Load = deltaP
R = DeltaP;

% Initial Norm
normR = norm(R);

% Initial Stiffness
if GRAD == 0
    K = K0;
elseif GRAD == 1
    K = Assemble_Secant_Stiffness(Q,input,type,ndsE,dofN);
elseif GRAD == 2
    K = Assemble_Secant_Stiffness(Q,input,type,ndsE,dofN);
    if n > 2
        K = Approx_Tangent(Pest(:,n-1),Pest(:,n-2),Q,Qsol(:,n-2),L,dofT);       % not working yet...
    end
end
fprintf('\n');
Kp = K; % Store previous stiffess matrix
% Display infor regarding stiffness matrix
fullK = full(K); 
maxK = max(max(K)); 
fprintf('\n Maximum stiffness value:\n'); disp(maxK);
minK = min(min(K));
fprintf('\n Minimum stiffness value:\n'); disp(minK);
diagK = diag(fullK);
fprintf('\n Diagonal terms:\n'); disp(diagK);
maxDiag = max(diagK);
fprintf('\n Maximum diagonal value:\n'); disp(maxDiag);
minDiag = min(diagK);
fprintf('\n Minimum diagonal value:\n'); disp(minDiag);
detK = det(K);
fprintf('\n Matrix Determinant:\n'); disp(detK);
detKboundary = det(K(L,L));
fprintf('\n Matrix determinant (with BCs):\n'); disp(detKboundary);

figure
spy(K)
figure
heatmap(fullK)


iter = 0;   % Count iterations

% Initialise Process Modification Arrays
IStage = 1;
CurrentStage = 1;
XO = zeros(dofT,1);
XN = zeros(dofT,1);
PhiN = zeros(dofT,1);
PhiO = zeros(dofT,1);
DiffN = zeros(dofT,1);
DiffO = zeros(dofT,1);
Beta = zeros(maxMOD,1);
LambdaO = 0; 
LambdaN = 0;
nNode = size(input.ND,1);

% Commence iterative process
while normR > alpha1

    % Begin line to print diagnostic data
    fprintf('\n\n   Iteration %i',iter);

    % Change in displacement estimate = K^-1 * R
    DeltaQ = SolveQ(K,R,L,input);
    fprintf('\n   Norm of DeltaQ = %.2d',norm(DeltaQ));

    % Update estimate
    Q = Q + DeltaQ';

    % If Process Modification is enabled...
    if PROCMOD == 1 && maxMOD > 1

        % Calculate an improved estimate
        Q = ImprovedGuess(Phi,Beta,Q,XO,IStage,dofT);
        
        % Assess Convergence, Modify Process
        [Q,XO,DiffN,DiffO,Phi,Beta,IStage,LambdaN,PhiN,...
        PhiO] = ProcessMod_Mario(Q,XO,DiffN,DiffO,Phi,Beta,dofT,IStage,...
        LambdaN,PhiN,PhiO,alpha2,maxMOD);
        fprintf('\n   Process Modification Data:')
        fprintf('\n      Lambda = %.5f',LambdaN);
        fprintf('\n      normPhi = %.2d',norm(PhiN));
        
    end 

    % Determine Load Vector at Q
    P = SolveP(Q,input,ndsE);
    fprintf('\n   Calculated load vector:\n');disp(P);
    
    % Calculate residual R
    R = full(Pn) - P';
    fprintf('\n   Calculated residual:\n');disp(R');
    
    % Normalise the residual vector
%         R = R/norm(R);

    % Calculate norm
%         normR = norm(R(L));

    % Mario's norm
    normR = sqrt(R * R')/(dofT*abs(max(Pn)));
    fprintf('\n   normResidual = %.2d',normR);

    % Reset Process Modification Arrays
    if IStage == CurrentStage+1
        CurrentStage = CurrentStage+1;
        PhiN = zeros(dofT,1);        % Process Modification Parameter 
        PhiO = zeros(dofT,1);        % Process Modification Parameter 
        DiffN = zeros(dofT,1);       % Process Modification Parameter 
        DiffO = zeros(dofT,1);       % Process Modification Parameter 
        fprintf('\n      PROCESS MODIFIED \n');
    else
        XO = Q;
    end

    % Update iteration count
    iter = iter + 1;

    % End loop if convergence is not occurring
    if iter > maxITER 
        fprintf('\n      DID NOT CONVERGE (reached maximum number of iterations)');
        FAIL = 1;
        fprintf('\n');
        fprintf('\n   Exiting Process Prematurely...');
        fprintf('\n      Successful Load Steps: %i',n-1);
        fprintf('\n      Percentage Complete:   %.2f%%',(n-1)/N*100);
        fprintf('\n');
        break;
    end

    % End loop if displacements are exploding
    if norm(Q) > 10^4
        fprintf('\n      DID NOT CONVERGE (displacements -> infinity)');
        FAIL = 1;
        fprintf('\n');
        fprintf('\n   Exiting Program Prematurely...');
        fprintf('\n      Successful Load Steps: %i',n-1);
        fprintf('\n      Percentage Completed:  %.2f%%',(n-1)/N*100);
        fprintf('\n');
        break;
    end

end
