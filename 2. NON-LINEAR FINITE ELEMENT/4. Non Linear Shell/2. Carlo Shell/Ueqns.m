% This script is used to derive symbolic expressions for the fundamental
% kinematic equations formulating the non-linear shell element derived
% alongside Carlo Sansour

% Created:      20 November 2017
% Last Updated: 11 January 2018

%% PRELIMINARIES

clear
clc
close all

%% DEFINE THE ELEMENT

% t = element length (X1 direction)     % not used yet...
% b = element length (X2 direction)     % not used yet...
% h = element thickness (z direction)
syms t b h
assume(t > 0)
assume(b > 0)
assume(h > 0)

%% DEFINE SYMBOLIC VECTORS:

% Vectors defining top and bottom surfaces xA(X1,X2) and xB(X1,X2)
syms xA1 xA2 xA3 xB1 xB2 xB3
xA = [xA1; xA2; xA3];
xB = [xB1; xB2; xB3];

% Scalar parameter lambda(X1,X2) and its derivatives w.r.t X1 and X2
syms lam lam_1 lam_2

% Vectors defining derivatives of xA and xB with respect to X1
syms xA1_1 xA2_1 xA3_1 xB1_1 xB2_1 xB3_1 
xA_1 = [xA1_1; xA2_1; xA3_1];
xB_1 = [xB1_1; xB2_1; xB3_1];

% Vectors defining derivatives of xA and xB with respect to X2
syms xA1_2 xA2_2 xA3_2 xB1_2 xB2_2 xB3_2
xA_2 = [xA1_2; xA2_2; xA3_2];
xB_2 = [xB1_2; xB2_2; xB3_2];

% Define rotation vector (three angles)
syms gam1 gam2 gam3
assume(gam1, 'real');
assume(gam2, 'real');
assume(gam3, 'real');
gam = [gam1; gam2; gam3];

%% DEPENDENT PARAMETERS

% Right stretch tensor U is a function of the following 24 variables:
var = [xA(1) xA(2) xA(3)...
       xB(1) xB(2) xB(3)...
       xA_1(1) xA_1(2) xA_1(3)...
       xB_1(1) xB_1(2) xB_1(3)...
       xA_2(1) xA_2(2) xA_2(3)...
       xB_2(1) xB_2(2) xB_2(3)...
       lam lam_1 lam_2...
       gam(1) gam(2) gam(3)];
   
%% CONSTANT PART OF DEFORMATION GRADIENT TENSOR   (i.e. F0)

% Vectors forming columns of constant part of Deformation Gradient (F0)

% "x0" = constant part of x

% dx0/dX1
x01_1 = symfun((1/2)*(-lam_1*h/2)*xA(1) + (1/2)*(1-lam*h/2)*xA_1(1) + (1/2)*(lam_1*h/2)*xB(1) + (1/2)*(1+lam*h/2)*xB_1(1), var);
x02_1 = symfun((1/2)*(-lam_1*h/2)*xA(2) + (1/2)*(1-lam*h/2)*xA_1(2) + (1/2)*(lam_1*h/2)*xB(2) + (1/2)*(1+lam*h/2)*xB_1(2), var);
x03_1 = symfun((1/2)*(-lam_1*h/2)*xA(3) + (1/2)*(1-lam*h/2)*xA_1(3) + (1/2)*(lam_1*h/2)*xB(3) + (1/2)*(1+lam*h/2)*xB_1(3), var);
x0_1 = {x01_1; x02_1; x03_1};

% dx0/dX2
x01_2 = symfun((1/2)*(-lam_2*h/2)*xA(1) + (1/2)*(1-lam*h/2)*xA_2(1) + (1/2)*(lam_2*h/2)*xB(1) + (1/2)*(1+lam*h/2)*xB_2(1), var);
x02_2 = symfun((1/2)*(-lam_2*h/2)*xA(2) + (1/2)*(1-lam*h/2)*xA_2(2) + (1/2)*(lam_2*h/2)*xB(2) + (1/2)*(1+lam*h/2)*xB_2(2), var);
x03_2 = symfun((1/2)*(-lam_2*h/2)*xA(3) + (1/2)*(1-lam*h/2)*xA_2(3) + (1/2)*(lam_2*h/2)*xB(3) + (1/2)*(1+lam*h/2)*xB_2(3), var);
x0_2 = {x01_2; x02_2; x03_2};

% dx0/dz
x01_3 = symfun((1/h)*xA(1) + (-1/h)*xB(1), var);
x02_3 = symfun((1/h)*xA(2) + (-1/h)*xB(2), var);
x03_3 = symfun((1/h)*xA(3) + (-1/h)*xB(3), var);
x0_3 = {x01_3; x02_3; x03_3};

F0 = {x0_1 x0_2 x0_3};

% Notation used to store functions in F0:

% F0{i}{j} == dx0i/dxj (i = row of F0, j = column of F0)

%% LINEAR PART OF DEFORMATION GRADIENT TENSOR   (i.e. F1)

% Vectors forming columns of linear part of Deformation Gradient (F1)

% "x1" = linear part of x

% dx1/dX1
x11_1 = symfun((1/h)*xA_1(1) + (-1/h)*xB_1(1), var);
x12_1 = symfun((1/h)*xA_1(2) + (-1/h)*xB_1(2), var);
x13_1 = symfun((1/h)*xA_1(3) + (-1/h)*xB_1(3), var);
x1_1 = {x11_1; x12_1; x13_1};

% dx1/dX2
x11_2 = symfun((1/h)*xA_2(1) + (-1/h)*xB_2(1), var);
x12_2 = symfun((1/h)*xA_2(2) + (-1/h)*xB_2(2), var);
x13_2 = symfun((1/h)*xA_2(3) + (-1/h)*xB_2(3), var);
x1_2 = {x11_2; x12_2; x13_2};

% dx1/dX3
x11_3 = symfun((2*lam/h)*xA(1) + (-2*lam/h)*xB(1), var);
x12_3 = symfun((2*lam/h)*xA(2) + (-2*lam/h)*xB(2), var);
x13_3 = symfun((2*lam/h)*xA(3) + (-2*lam/h)*xB(3), var);
x1_3 = {x11_3; x12_3; x13_3};

F1 = {x1_1 x1_2 x1_3};

% Notation used to store functions in F1:

% F1{i}{j} == dx1i/dxj (i = row of F1, j = column of F1)

%% CONSTANT PART OF RIGHT STRETCH TENSOR   (i.e. U0)

U01 = cell(3,1);    % pre-allocate size
U02 = cell(3,1);
U03 = cell(3,1);
delta = [1 0 0; 0 1 0; 0 0 1]; % Kronecker Delta
for i = 1:3
    for j = 1:3
        xTerm1 = F0{i}{j};
        
        xTerm2 = 0;
        for k = 1:3
            rTerm2 = (1-cos(norm(gam)))/(norm(gam)^2)*(gam(i)*gam(k)-dot(gam,gam)*delta(i,k));
            
            rTerm1 = 0;
            for s = 1:3
                rTerm1 = rTerm1 + -sin(norm(gam))/norm(gam)*eijk([i;j;s],1)*gam(s);
            end
            xTerm2 = xTerm2 + F0{k}{j}*(rTerm1+rTerm2);
            
        end
        if i == 1;  U01{j} = xTerm1 + xTerm2; end
        if i == 2;  U02{j} = xTerm1 + xTerm2; end
        if i == 3;  U03{j} = xTerm1 + xTerm2; end
        
    end
end
U0 = {U01 U02 U03};

%% LINEAR PART OF RIGHT STRETCH TENSOR U   (i.e. K)

U11 = cell(3,1);    % pre-allocate size
U12 = cell(3,1);
U13 = cell(3,1);
delta = [1 0 0; 0 1 0; 0 0 1]; % Kronecker Delta
for i = 1:3
    for j = 1:3
        xTerm1 = F1{i}{j};
        
        xTerm2 = 0;
        for k = 1:3
            rTerm2 = (1-cos(norm(gam)))/(norm(gam)^2)*(gam(i)*gam(k)-dot(gam,gam)*delta(i,k));
            
            rTerm1 = 0;
            for s = 1:3
                rTerm1 = rTerm1 + -sin(norm(gam))/norm(gam)*eijk([i;j;s],1)*gam(s);
            end
            xTerm2 = xTerm2 + F1{k}{j}*(rTerm1+rTerm2);
            
        end
        if i == 1;  U11{j} = xTerm1 + xTerm2; end
        if i == 2;  U12{j} = xTerm1 + xTerm2; end
        if i == 3;  U13{j} = xTerm1 + xTerm2; end
        
    end
end
K = {U11 U12 U13};

% This loop can be combined with the constant part loop to save time. It
% has been left like this for now just for clarity of the distinct steps
% involved in the formulation of this element

% Constant and Linear parts form U = U0 + z*K

%% FIRST VARIATION OF RIGHT STRETCH TENSOR U

% Requires the differentiation of U0 and K with respect to all of their 
% dependent variables. There are a total of 24 dependent parameters for
% each tensor

% U0 and K consist of 9 components each. Differentiation of each
% component with respect each of the 24 parameters yields a total of:
% 216 = 9 x 24 expressions per tensor (majority of which will be zero)

% Altogether there are 432 = 2 x 216 expressions to store. It was found
% that of these, 198 are non-zero expressions.

% Partial derivatives of components are able to be found using:

%   dU_ij/dm = diff(U{i}{j},m)

num_param = numel(var);     % number of dependent parameters

dU0 = cell(3,3,num_param);  % pre-allocate size
dK  = cell(3,3,num_param);  % pre-allocate size
nu = 0;
nk = 0;
for a = 1:num_param
    par = var(a);
    
    for i = 1:3
        for j = 1:3
            
        % Differentitate component of U0 w.r.t current scalar parameter
            dU0_p = diff(U0{i}{j},par);
            
        % Differentitate component of K  w.r.t current scalar parameter
            dK_p  = diff(K{i}{j}, par);
            
        % Store expressions
            dU0{i}{j}{a} = dU0_p;
             dK{i}{j}{a} = dK_p;
        
        % Count non-zero entries:
            if dU0_p == 0; nu = nu + 1; end
            if dK_p  == 0; nk = nk + 1; end
                
        end
    end
end

fprintf('\n   Number of non-zero expression in U0: %i',216-nu);
fprintf('\n   Number of non-zero expression in K:  %i',216-nk);
% Note that expressions are stored as:

% dU0_ij/dm = dU0{i)(j){m)
%  dK_ij/dm =  dK{i}{j}{m}


%% SECOND VARIATION OF RIGHT STRETCH TENSOR U

% Requires differentiation of U0 and K to all unique "pairs" or parameters.
% There are 276 unique pairs of parameters. Therefore, differentiting the 9
% components within either U0 and K to each ofthe 276 pairs yields 2484
% expressions. Altogether we have 2 x 2484 = 4968 expressions, the majority
% of which will be zero.

% If the order of differentiation of each pair of parameters is significant
% then there are 552 second order derivatives to find for each component
% giving a total of 9936 = 552 x 9 x 2

% Mixed derivatives can be found using:

%   dU_ij/dm1dm2 = diff(U{i}{j},m1,m2)

ddU0 = cell(3,3,num_param,num_param);
ddK  = cell(3,3,num_param,num_param);
nu = 0; nk = 0;
for a = 1:num_param
    par1 = var(a);
    
    for b = 1:num_param
        par2 = var(b);
        
        if a >= b
            dU0_p1p2 = diff(U0{i}{j},par1,par2);
            dK_p1p2  = diff(K{i}{j} ,par1,par2);
        
            ddU0{i}{j}{a}{b} = dU0_p1p2;
             ddK{i}{j}{a}{b} = dK_p1p2;
             
            if dU0_p1p2 == 0; nu = nu + 1; end
            if dK_p1p2  == 0; nk = nk + 1; end

        end
                       
    end
end

fprintf('\n   Number of non-zero expression in U0: %i',2484-nu);
fprintf('\n   Number of non-zero expression in K:  %i',2484-nk);

% Once understood, this can be inserted within the 'first derivatives' loop
% as an extra section. Could be made more efficient by only performing
% differentiation if first derivative is non-zero (as all mixed
% derivatives involving that variable will also be zero)