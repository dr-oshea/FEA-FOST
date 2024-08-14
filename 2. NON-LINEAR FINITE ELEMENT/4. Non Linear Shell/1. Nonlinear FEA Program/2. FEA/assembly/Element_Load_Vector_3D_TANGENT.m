function Pe = Element_Load_Vector_3D_TANGENT (q,elCoords,theta,model_data,fibres)

% This script generates the element load vector at the displacement field
% given by Q for a 3D brick element. This is a nonlinear relationship. The 
% integral is determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  04 April 2023

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
I = eye(3);

% Global Element Constitutive Law
[~,~,T,~,~] = Transform2(phi);
% Lambda = inv(T) * Lambda * T ;
% Mu     = inv(T) *   Mu   * T ;


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
            
            % Displacement Gradient Tensor
            GradU = V2M(B * q);

            % Deformation gradient vector
            F = GradU + I;
            
            % Assemble matrix for fourth-order symmetry projection tensor:
            S4 = (1/2)*(circledot22(I,I) + circlestar22(I,I));
            
            % Assemble (I odot GradU) tensor:
            A1 = DblCon44(S4,circledot22(F.',I));
            A1m = T2M(A1);

            % Get PK2 stress tensor, and tangent stiffness tensor           
            [S,~] = Get_S_CC(F, model_data, fibres);
        
            % Load vector at a point
            P = B.' * A1m.' * M2V(S) .* det(J);
        
            % Update Gaussian sum (integral)
            Pe = Pe + (int_weight(i) * int_weight(j) * int_weight(k)) .* P;
            
        end

    end
   
end

