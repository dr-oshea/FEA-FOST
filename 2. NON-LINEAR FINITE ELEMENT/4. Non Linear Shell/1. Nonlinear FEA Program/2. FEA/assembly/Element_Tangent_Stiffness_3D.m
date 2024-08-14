function Kele = Element_Tangent_Stiffness_3D (q, elCoords, angleDegrees, modelData, fibres)

% This script generates the element tangent stiffness matrix at a 
% given displacement field for a 3D brick element. The integral is 
% determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  23 March 2023

% INPUTS:
% q = local vector of known displacements
% elCoords = nodal cordinates for local nodes
% theta = material orientation

% OUTPUTS:
% K0e = element stiffness matrix

%% ---------------------------------------------------------------------------
nodesPerElement = 8;
dofPerNode = 3;
angleRadians = deg2rad(angleDegrees);
I = Iden2(3);

% Global Element Constitutive Law
[~,~,T9,~,~] = Transform2(angleRadians);

%% GAUSSIAN INTEGRATION

% Retrieve Gauss points (integration weights and points)
[intPoint,intWeight] = GaussPoints(2);

% Initialise element stiffness matrix K
Kele  = zeros(nodesPerElement*dofPerNode,nodesPerElement*dofPerNode); 

for i = 1:numel(intPoint)
     
    xi = intPoint(i);
   
    for j = 1:numel(intPoint)
        
        eta = intPoint(j);
        
        for k = 1:numel(intPoint)
        
            mu = intPoint(k);
            
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
            
            % Assemble (I odot GradU) tensor, convert to matrix
            A1 = DblCon44(S4,circledot22(F.',I));
            A1m = T2M(A1);

            % Get PK2 stress tensor, and tangent stiffness tensor           
            [S,CC] = Get_S_CC(F, modelData, fibres);
            
            Dm = T2M(circledot22(I,S)) + A1m.'*T2M(CC)*A1m;
            
            % Calcualte element tangent stiffness in global coords
            K = (B' * Dm * B) .* det(J);

            Kele = Kele + (intWeight(i) * intWeight(j) * intWeight(k)) .* K ;
        
        end
    end
end
