function Pe = Element_Load_Vector_3D (q,elCoords,theta,Lambda,Mu)

% This script generates the element load vector at the displacement field
% given by Q for a 3D brick element. This is a nonlinear relationship. The 
% integral is determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  13 April 2018

% INPUTS:
% q = local displacement vector
% elCoords = nodal cordinates for local nodes
% theta = material orientation
% Lambda = Lambda material tensor (flattened) - material coordinates
% Mu = Mu material tensor (flattened) - material coordinates

% OUTPUTS:
% Pe = element load vector

%% ---------------------------------------------------------------------------

ndsPerElem = 8;
dofN = 3;
phi = deg2rad(theta);

% Global Element Constitutive Law
[~,~,T,~,~] = Transform2(phi);
Lambda = inv(T) * Lambda * T ;
Mu     = inv(T) *   Mu   * T ;


%% Gaussian Integration

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(2);

% Initialise K and B matrices
Pe  = zeros(ndsPerElem*dofN,1); 
% B   = zeros(dofN*dofN,ndsPerElem*dofN);

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

            % Compute stress vector:
            Svec = S_Vector(Cvec,Lambda,Mu);
        
            % Load at a point
            P = (1/2)* transpose(B) * transpose(A1 + A3) * Svec .* det(J);
        
            % Update Gaussian sum (integral)
            Pe = Pe + (int_weight(i) * int_weight(j) * int_weight(k)) .* P;
            
        end

    end
   
end

