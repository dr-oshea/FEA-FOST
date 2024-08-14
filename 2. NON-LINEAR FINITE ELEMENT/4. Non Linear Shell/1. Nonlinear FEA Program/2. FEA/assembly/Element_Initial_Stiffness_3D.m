function K0e = Element_Initial_Stiffness_3D (elCoords, theta, model_data, fibres)

% This script generates the element stiffness matrix in the linear region
% for a 3D brick element. The integral is determined using Gaussian
% integration

% Author:   Daniel O'Shea
% Created:  12 April 2018

% INPUTS:
% elCoords = nodal cordinates for local nodes
% theta = material orientation
% C = infinitesimal material tensor (flattened) in material coordinates

% OUTPUTS:
% K0e = element stiffness matrix

%% ---------------------------------------------------------------------------

ndsPerElem = 8;
dofN = 3;

phi = deg2rad(theta);

% Global Element Constitutive Law
[~,C] = Get_S_CC(eye(3),model_data, fibres);
C = T2M(C);

[~,~,T,~,~]=Transform2(phi);
Cg = inv(T) * C * T ;

%% GAUSSIAN INTEGRATION

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(3);

% Initialise K and B matrices
K0e  = zeros(ndsPerElem*dofN,ndsPerElem*dofN); 

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

            A1 = [2 0 0 0 0 0 0 0 0;
                  0 2 0 0 0 0 0 0 0;
                  0 0 2 0 0 0 0 0 0;
                  0 0 0 1 0 0 1 0 0;
                  0 0 0 0 1 0 0 1 0;
                  0 0 0 0 0 1 0 0 1;
                  0 0 0 1 0 0 1 0 0;
                  0 0 0 0 1 0 0 1 0;
                  0 0 0 0 0 1 0 0 1];
                  
            % Strain Displacement Relationship
            B = B_Matrix_3D(dN);

            B0 = (1/2).*(A1 * B);

            K = (B0' * Cg * B0) .* det(J);

            K0e = K0e + int_weight(i)*int_weight(j)*int_weight(k).*K;
            
        end

    end
   
end

% point_niu = [ -1  -1   1   1]; 
% point_xi  = [ -1   1   1  -1]; 
% 
% 
% % Store B1-matrices at nodes (for post processing)
% for i=1:numel(point_niu)
%     % Shape functions (Natural Coords):
%     N1 = 1/4*(1-xi)*(1-niu);
%     N2 = 1/4*(1+xi)*(1-niu);
%     N3 = 1/4*(1+xi)*(1+niu);
%     N4 = 1/4*(1-xi)*(1+niu);
%     Nmat = [N1; N2; N3; N4];
% 
%     % Shape function derivatives (Natural Coords):
%     R = 1/4.*[(-1+niu) (1-niu) (1+niu) (-1-niu); 
%                (-1+xi) (-1-xi)  (1+xi)   (1-xi)];
% 
%     % Jacobian Matrix
%     J =  R * elCoords;
% 
%     % Shape function derivatives (Global Coords):
%     dN = inv(J)*R;
% 
%     A1 = [2 0 0 0;
%           0 2 0 0;
%           0 0 1 1;
%           0 0 1 1];
% 
%     % Strain Displacement Relationship
%     B = B1_Matrix_2D(dN);
% 
%     B(1:dofN*dofN,:,i+1) = B; 
%     
% end
