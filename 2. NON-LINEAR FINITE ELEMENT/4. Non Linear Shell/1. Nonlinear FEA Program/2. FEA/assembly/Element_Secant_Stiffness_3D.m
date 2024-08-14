function Ke = Element_Secant_Stiffness_3D (q, elCoords, theta, L, M)

% This script generates the element secant stiffness matrix in the at a 
% given displacement field for a 3D brick element. The integral is 
% determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  12 April 2018

% INPUTS:
% q = local vector of known displacements
% elCoords = nodal cordinates for local nodes
% theta = material orientation
% C = infinitesimal material tensor (flattened) in material coordinates

% OUTPUTS:
% Ke = element stiffness matrix

%% ---------------------------------------------------------------------------
global nset cset
ndsPerElem = 8;
dofN = 3;
phi = deg2rad(theta);
I = Iden2(3);
num_param = numel(nset);

% Global Element Constitutive Law
[~,~,T,~,~] = Transform2(phi);
Lg = inv(T) * L * T ;
Mg = inv(T) * M * T ;

Lg = M2T(Lg,'dot');
Mg = M2T(Mg,'dot');

%% GAUSSIAN INTEGRATION

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(3);

% Initialise K and B matrices
Ke  = zeros(ndsPerElem*dofN,ndsPerElem*dofN); 

for i = 1:numel(int_point)
     
    xi = int_point(i);
   
    for j = 1:numel(int_point)
        
        eta = int_point(j);
        
        for k = 1:numel(int_point)
            
            mu = int_point(k);
        
            % Shape function derivatives (Natural Coords):
            R = 1/8.*[-(1-eta)*(1-mu) (1-eta)*(1-mu) (1+eta)*(1-mu) -(1+eta)*(1-mu) -(1-eta)*(1+mu) (1-eta)*(1+mu) (1+eta)*(1+mu) -(1+eta)*(1+mu);...
               -(1-xi)*(1-mu) -(1+xi)*(1-mu) (1+xi)*(1-mu) (1-xi)*(1-mu) -(1-xi)*(1+mu) -(1+xi)*(1+mu) (1+xi)*(1+mu) (1-xi)*(1+mu);...
               -(1-xi)*(1-eta) -(1+xi)*(1-eta) -(1+xi)*(1+eta) -(1-xi)*(1+eta) (1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

            % Jacobian Matrix
            J =  R * elCoords;

            % Shape function derivatives (Global Coords):
            dN = inv(J)*R;

            % Displacement Gradient Relationship
            B = B_Matrix_3D(dN);

            % Displacement Gradient Vector
            GradU = B * q;

            % Assemble A1 matrix:
            A1 = [2 0 0 0 0 0 0 0 0;
                  0 2 0 0 0 0 0 0 0;
                  0 0 2 0 0 0 0 0 0;
                  0 0 0 1 0 0 1 0 0;
                  0 0 0 0 1 0 0 1 0;
                  0 0 0 0 0 1 0 0 1;
                  0 0 0 1 0 0 1 0 0;
                  0 0 0 0 1 0 0 1 0;
                  0 0 0 0 0 1 0 0 1];

            % Assemble A2 matrix:
            A2 = [GradU(1)     0        0        0    GradU(5)     0        0        0    GradU(9) ;
                      0    GradU(2)     0        0        0    GradU(6) GradU(7)     0        0    ;
                      0        0    GradU(3) GradU(4)     0        0        0    GradU(8)     0    ;
                      0        0    GradU(7) GradU(2)     0        0        0    GradU(6)     0    ;
                  GradU(8)     0        0        0    GradU(3)     0        0        0    GradU(4) ;
                      0    GradU(9)     0        0        0    GradU(1) GradU(5)     0        0    ;
                      0    GradU(4)     0        0        0    GradU(8) GradU(3)     0        0    ;
                      0        0    GradU(5) GradU(9)     0        0        0    GradU(1)     0    ;
                  GradU(6)     0        0        0    GradU(7)     0        0        0    GradU(2) ];

            % Assemble A3 matrix:
            A3=[2*GradU(1)     0         0          0    2*GradU(5)     0          0         0    2*GradU(9);
                     0    2*GradU(2)     0          0         0    2*GradU(6) 2*GradU(7)     0         0    ;
                     0         0    2*GradU(3) 2*GradU(4)     0         0          0    2*GradU(8)     0    ;
                     0      GradU(4)  GradU(7)   GradU(2)     0      GradU(8)   GradU(3)  GradU(6)     0    ;
                  GradU(8)     0      GradU(5)   GradU(9)  GradU(3)     0          0      GradU(1)  GradU(4);
                  GradU(6)  GradU(9)     0          0      GradU(7)  GradU(1)   GradU(5)     0      GradU(2);
                     0      GradU(4)  GradU(7)   GradU(2)     0      GradU(8)   GradU(3)  GradU(6)     0    ;
                  GradU(8)     0      GradU(5)   GradU(9)  GradU(3)     0          0      GradU(1)  GradU(4);
                  GradU(6)  GradU(9)     0          0      GradU(7)  GradU(1)   GradU(5)     0      GradU(2)];

            % Right Cauchy-Green deformation vector:
            Cvec = [1;1;1;0;0;0;0;0;0] + (A1 + A2)*GradU;

            % Right Cauchy-Green deformation tensor:
            Cmat = V2M(Cvec);

            % Simplified form of strain derivatives (small displacement)
            CC_lam = DblCon44(circledot22(Power2(Cmat,-1),I),Lg,'dot');
            CC_mu = zeros(3,3,3,3);
            for m = 1:num_param
                ni = nset(m);
                ci = cset(m);
                CC_mu = CC_mu + DblCon44(circledot22(Power2(Cmat,ni/2-1),I),ci*Mg,'dot');
            end
            Cg = T2M(CC_lam + CC_mu, 'dot');

            % Element secant stiffness
            K = (1/4) .* (B' * (A1+A3)' * Cg * (A1+A2) * B) .* det(J) ;

            Ke = Ke + (int_weight(i) * int_weight(j) * int_weight(k)) .* K ;
            
        end

    end
   
end
