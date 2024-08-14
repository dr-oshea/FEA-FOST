function Pe = Element_Load_Vector_2D_TANGENT (q,elCoords,theta,model_data)

% This script generates the element load vector at the displacement field
% given by Q for a 2D quad element. This is a nonlinear relationship. The 
% integral is determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  19 June 2023

% INPUTS:
% q = local displacement vector
% elCoords = nodal cordinates for local nodes
% theta = material orientation
% Lambda = Lambda material tensor (flattened) - material coordinates
% Mu = Mu material tensor (flattened) - material coordinates

% OUTPUTS:
% Pe = element load vector

%% ---------------------------------------------------------------------------
ndsPerElem = 4;
dofN = 2;
phi = deg2rad(theta);
I = eye(2);

% Global Element Constitutive Law
[T4,~,~,~,~] = Transform2(phi);
% Lambda = inv(T4) * Lambda * T4 ;
% Mu     = inv(T4) *   Mu   * T4 ;


%% Gaussian Integration

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(3);

% Initialise K and B matrices
Pe  = zeros(ndsPerElem*dofN,1); 
% B   = zeros(dofN*dofN,ndsPerElem*dofN);

for i = 1:numel(int_point)
     
    xi = int_point(i);
   
    for j = 1:numel(int_point)
        
        eta = int_point(j);
        
        % Shape function derivatives (Natural Coords):
        R = 1/4.*[(-1+eta) (1-eta) (1+eta) (-1-eta); 
                   (-1+xi) (-1-xi)  (1+xi)   (1-xi)];
        
        % Jacobian Matrix
        J =  R * elCoords;
        
        % Shape function derivatives (Global Coords):
        dN = inv(J)*R;
        
        % Displacement Gradient Relationship
        B = B_Matrix_2D(dN);

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
        [S,~] = Get_S_CC(F, model_data);
    
        % Load vector at a point
        P = B.' * A1m.' * M2V(S) .* det(J);
        
        % Update Gaussian sum (integral)
        Pe = Pe + int_weight(i) * int_weight(j) .* P;

    end
   
end

