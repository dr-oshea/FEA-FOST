%% Testing Process Modification Program on 1D problem

clear
close all
clc

%% Define the the function to model f(x)
syms x

% Taylor series of log function
% f(x) = taylor(log(x),x,1,'Order',4)+ (1/1+1/2+1/3);

% Polynomial function
a = 2; b = 1;
f(x) = a*x^2+b*x;

% Exponential function
% f(x) = exp(x)-1;

%% Plot the known function
xvec = 0:0.01:1;
yvec = f(xvec);
figure
hold on
plot(xvec,yvec)

%% Derivative of function f'(x)
df(x) = diff(f(x),x);

%% Initialise "load" steps for function

% Inital and final function values
P0 = 0;
Pf = 2;

% Define total divisions of "load"
N = 20;

% Calculate "load" step
DeltaP = (Pf-P0)/(N);

% Set of "load steps to run model at
f_vec = P0:DeltaP:Pf;
numf = numel(f_vec);

%% Initial Gradient Method
GRAD = 1;
PROCMOD = 1;


% Initial Gradient (derivaive of function at initial x value)
Q0 = 0;
initial_grad = double(df(Q0));

% Define tolerance for residual
tol = 10^-4;

dofT = size(Q0,1);
Qsol = zeros(dofT,N);
Psol = zeros(dofT,N);
maxMOD = 1;

for n = 1:N
    
    fprintf('\n\n   COMMENCING LOAD STEP #%i',n);
    if n == 1
        Q = Q0;
        P = P0;
    else
        Q = Qsol(n-1);
        P = Psol(n-1);
    end
    
    Pn = P + DeltaP;
    
    R = Pn - P;
    
    normR = norm(R);
    
    if GRAD == 0
        K = initial_grad;
    elseif GRAD == 1
        if n == 1
            K = initial_grad;
        else
            K = (P-P0)/(Q-Q0);
        end
    end
    
    % Initialise Variables for Iterative Process
    iter = 0;
    subiter = 1;
    mod = 0;
    Beta = zeros(maxMOD,1);
    Phi0 = zeros(dofT,1);
    Lambda0 = 0;
    DeltaQ0 = zeros(dofT,1);
    alpha2 = 10^-2;
    
    fprintf('\n      Iter  |   norm(DelQ)  |   max(Q)   |   Lambda   |  norm(Phi)  |  norm(R)');
    
    while normR > tol
        
        fprintf('\n      %i',iter);
        
        % Change in guess
        DeltaQ = K^-1 * R;
        fprintf('         %.2d',norm(DeltaQ));
        
        % Updated guess
        Q = Q + DeltaQ;
        
        if PROCMOD == 1 && mod >= 1 && mod <= dofT
            
            % Current guess
            Current = Q;
            
            % Modification to guess
            Term = zeros(dofT,1);
            Modified = zeros(dofT,1);
            for i = 1:mod
                Term = Beta(i)*Phi*(DeltaQ+Phi*Term');
                Modified = Modified + Term;
            end
            
            % Improved guess
            Q = Current + Modified;
            
        end
        fprintf('       %.2d',max(Q));
        
        % Calculate function and residual
        P = double(f(Q));
        
        R = Pn - P;
        
        normR = norm(R);
        
        % Potentially modify the iterative process
        if subiter >= 2
            
            Lambda = (DeltaQ*DeltaQ0')/(DeltaQ0*DeltaQ0');
            Phi = DeltaQ0/sqrt(DeltaQ0*DeltaQ0');
            
            if Lambda < 0; fprintf('     %.5f',Lambda); end
            if Lambda >= 0; fprintf('      %.5f',Lambda); end
            fprintf('     %7.2e',norm(Phi));
            
            if (dofT == 1 && mod < maxMOD) || (norm(Lambda-Lambda0)<alpha2 && norm(Phi-Phi0)<alpha2 && mod < maxMOD)
                Beta(mod+1) = Lambda/(1-Lambda);
                mod = mod + 1;
                subiter = 0;
            else
                Lambda0 = Lambda;
                Phi0 = Phi; 
            end
            
        else
            fprintf('      -------      -------');

        end
        
        fprintf('      %.2d',normR);

        if subiter == 0; fprintf('\n      PROCESS MODIFIED'); ON = 0; end
        DeltaQ0 = DeltaQ;
        
        iter = iter + 1;
        subiter = subiter + 1;
        
    end
    
    Qsol(n) = Q;
    Psol(n) = P;
    
    
end
fprintf('\n');

scatter(Qsol,Psol)