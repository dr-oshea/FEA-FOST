%% ABOUT THIS SCRIPT

% This script constucts a non-linear shell element, matching the
% formulation provided in "Ueqns.m"

% Created:          November 2017
% Last Updater:     November 2017


% It is a work in progress that needs to be updated....

%% PRELIMINARIES

clear
clc
close all

%% DEFINE THE ELEMENT

% t = xdim
% b = ydim
% h = zdim
syms t b h
assume(t > 0)
assume(b > 0)
assume(h > 0)

% Nodal Coordinates (Reference Configuration)
x1 = -t/2; y1 = -b/2; z1 = -h/2;
x2 = +t/2; y2 = -b/2; z2 = -h/2;
x3 = +t/2; y3 = +b/2; z3 = -h/2;
x4 = -t/2; y4 = +b/2; z4 = -h/2;
x5 = -t/2; y5 = -b/2; z5 = +h/2;
x6 = +t/2; y6 = -b/2; z6 = +h/2;
x7 = +t/2; y7 = +b/2; z7 = +h/2;
x8 = -t/2; y8 = +b/2; z8 = +h/2;
xbot = [x1;x2;x3;x4];
xtop = [x5;x6;x7;x8];
ybot = [y1;y2;y3;y4];
ytop = [y5;y6;y7;y8];


%% FOR A GIVEN POINT ON THE X-Y PLANE (at any z... -> we integrate over z later...)
syms xi eta

% Shape functions (2D plane)
N1 = 1/4*(1-xi)*(1-eta);  
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);
Nmat = [N1 N2 N3 N4];

% Derivatives of Shape Functions w.r.t Natural Coordinates
R = 1/4.*[(-1+eta) (1-eta) (1+eta) (-1-eta); 
          (-1+xi)  (-1-xi) (1+xi)  (1-xi)  ];  % [d.N/d.xi ; d.N/d.eta]

% Jacobian: Derivatives of Global Coordinates w.r.t Natural Coordinates
J =  R * [xtop ytop];     %  [d.x/d.xi , d.y/d.xi ; d.x/d.eta d.y/d.eta]

% Derivative of Shape Functions w.r.t Global Coordinates
dN = inv(J)*R;    % = [d.N/d.x; d.N/d.y];

% X1 and X2 (Position Vectors of top and bottom surface)
X1 = dot([Nmat;Nmat],[xtop ytop]',2);
X2 = dot([Nmat;Nmat],[xbot ybot]',2);

% Derivatives of X1 and X2
X1_1 = dot([dN(1,:);dN(1,:)],[xtop ytop]',2);
X1_2 = dot([dN(2,:);dN(2,:)],[xtop ytop]',2);
X2_1 = dot([dN(1,:);dN(1,:)],[xbot ybot]',2);
X2_2 = dot([dN(2,:);dN(2,:)],[xbot ybot]',2);

% Guess Lambda
sym lambda1 lambda2 lambda3 lambda4
lam = [lambda1 lambda2 lambda3 lambda4];
lambda = dot(Nmat,lam);
lambda_1 = dot(dN(1,:),lam);
lambda_2 = dot(dN(2,:),lam);

% Vectors Forming Columns of X0 part of Deformation Gradient
X0_1 = (1/2)*(-lambda_1*h/2)*X1 + (1/2)*(1-lambda*h/2)*X1_1 + (1/2)*(lambda_1*h/2)*X2 + (1/2)*(1+lambda*h/2)*X2_1;
X0_2 = (1/2)*(-lambda_1*h/2)*X1 + (1/2)*(1-lambda*h/2)*X1_2 + (1/2)*(lambda_1*h/2)*X2 + (1/2)*(1+lambda*h/2)*X2_2;
X0_3 = (1/h)*X1 + (-1/h)*X2;

        