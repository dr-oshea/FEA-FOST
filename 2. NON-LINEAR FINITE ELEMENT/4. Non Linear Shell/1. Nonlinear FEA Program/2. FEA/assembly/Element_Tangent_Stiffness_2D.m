function Ke = Element_Tangent_Stiffness_2D (q,elCoords, theta, model_data)

% This script generates the element tangent stiffness matrix at a 
% given displacement field for a 2D quad element. The integral is 
% determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  22 March 2022

% INPUTS:
% q = local vector of known displacements
% elCoords = nodal cordinates for local nodes
% theta = material orientation

% OUTPUTS:
% K0e = element stiffness matrix

%% ---------------------------------------------------------------------------
ndsPerElem = 4;
dofN = 2;
phi = deg2rad(theta);
I = Iden2(2);
IoI = circledot22(I,I);

% Global Element Constitutive Law
[T4,~,~,~,~] = Transform2(phi);


%% GAUSSIAN INTEGRATION

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(2);

% Initialise K and B matrices
Ke  = zeros(ndsPerElem*dofN,ndsPerElem*dofN); 

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
        [S,CC] = Get_S_CC(F, model_data);
        
        Dm = T2M(circledot22(I,S)) + A1m.'*T2M(CC)*A1m;
        
        % Calcualte element tangent stiffness in global coords
        K = (B' * Dm * B) .* det(J);
        
        Ke = Ke + (int_weight(i) * int_weight(j)) .* K ;

    end
    hello = 1;
    
    
end
