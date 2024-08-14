%% Process Modification - Testing on a various function

% This script tests the initial gradient method on various functions and
% applied process modification to improve the convergence

clc
close all
clear

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
fi = 0;
ff = 2;

% Define total divisions of "load"
div = 20;

% Calculate "load" step
delta_f = (ff-fi)/(div);

% Set of "load steps to run model at
f_vec = fi:delta_f:ff;
numf = numel(f_vec);

%% Initial Gradient Method

% Initial Gradient (derivaive of function at initial x value)
x0 = 0;
initial_grad = double(df(0));

% Define tolerance for residual
tol = 10^-4;

% Run model
xsol = zeros(numf,1);   % array to store solved values of x at each load step
count = zeros(numf,1);   % array to store no. of iterations at each load step
for n = 2:numf   % "for each load step..."
    
    % Previous load at known solution/iteration
    f0 = f_vec(n-1);
    
    % Current load at unknown solution/iteration
    f1 = f0 + delta_f;
    
    % Current residual
    R0 = f1 - f0;
    
    if n == 2       % for first load step, initial solution is initial x
        xs = x0;
    else            % for all other load steps, initial solution is previous known solution
        xs = xsol(n-1); 
    end
    
    % Initial guess;
    xi0 = xs;
    
    % First guess:
    dxi1 = initial_grad^-1 * R0;
    xi1 = xi0 + dxi1;
    
    % Second guess
    fi1 = double(f(xi1));
    R1 = f1 - fi1;
    dxi2 = initial_grad^-1 * R1;
    xi2 = xi1 + dxi2;
    
    % Process modification
    lambda = (xi2-xi1)/(xi1-xi0);
    beta = lambda/(1-lambda);
    
    xs = xi2 + beta*(xi2-xi1);
    R = f1 - double(f(xs));
    count(n-1) = 3;     % three guesses so far
    while abs(R) > tol
        
        dx = initial_grad^-1 * R;
        xs = xs + dx + beta * dx;
        
        R = f1 - double(f(xs));
        count(n-1) = count(n-1) + 1;
    end
    
    % Once tolerance is satisfied, store the solution of x at load step
    xsol(n) = xs;
    % Plot the point
    scatter(xs,f1)
    
end
hold off

% Plot number of iterations per load step
figure
plot(f_vec,count)