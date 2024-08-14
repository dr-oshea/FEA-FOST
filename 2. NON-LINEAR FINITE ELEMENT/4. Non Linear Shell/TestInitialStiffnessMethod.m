%% Initial Gradient Method - Testing on various scalar functions

% This script tests the initial gradient method on the logarithmic
% function. This will later be applied to the finite element method to
% determine a given force-displacement curve defined by minimum potential
% energy

clc
close all
clear

% Gradient Method: 1 = Initial, 2 = Secant, 3 = Tangent
METHOD = 2;

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
xlabel('x')
ylabel('f(x)')

%% Derivative of function f'(x)
df(x) = diff(f(x),x);

%% Initialise "load" steps for function

% Inital and final function values
f_i = 0;
f_f = 2;

% Define total divisions of "load"
div = 20;

% Calculate "load" step
delta_f = (f_f-f_i)/(div);

% Set of "load steps to run model at
f_vec = f_i+delta_f:delta_f:f_f;
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
for n = 1:numf   % "for each load step..."
    
    % Starting function value
    if n == 1
        f0 = f_i;
        xs = x0;
        k = initial_grad;
    else
        f0 = f_vec(n-1);
        xs = xsol(n-1);
        
        if METHOD == 1
            k = initial_grad;
        elseif METHOD == 2
            k = (f0-f_i)/(xs-x0);
            f_i = f0;
            x0 = xs;
        elseif METHOD == 3
            k = double(df(xsol(n-1)));
        end
    
    end
    
    % targeted function value
    f1 = f0 + delta_f;
    
    % Current residual
    R = f1 - f0;
    
    count(n) = 1;  % used to count number of iterations
    
    % "While the residual is larger than the defined tolerance..."
    while abs(R) > tol
        
        % Calculate change in x
        dx = k^-1 * R;
        
        % Update guess for x
        xi = xs + dx;
        
        % Determine function at guess
        fi = double(f(xi));
        
        % Calculate residual from goal load
        R = f1 - fi;
        
        % Update solution
        xs = xi;

        count(n) = count(n) + 1;

    end
    
    % Once tolerance is satisfied, store the solution of x at load step
    xsol(n) = xs;
    % Plot the point
    scatter(xs,f1)
    
end
hold off


% Plot number of iterations per load step
figure
plot(f_vec,count,'-x')
xlabel('f(x)')
ylabel('Number of Iterations')