function K0e = Element_Initial_Stiffness_2D (elCoords, theta, model_data)

% This script generates the element stiffness matrix in the linear region
% for a 2D quad element. The integral is determined using Gaussian
% integration

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% elCoords = nodal cordinates for local nodes
% theta = material orientation
% C = infinitesimal material tensor (flattened) in material coordinates

% OUTPUTS:
% K0e = element stiffness matrix

%% ---------------------------------------------------------------------------

ndsPerElem = 4;
dofN = 2;
phi = deg2rad(theta);

% Global Element Constitutive Law
[~,C] = Get_S_CC(zeros(2),model_data);
C = T2M(C);

[T4,~,~,~,~]=Transform2(phi);
Cg = inv(T4) * C * T4 ;

%% GAUSSIAN INTEGRATION

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(3);

% Initialise K and B matrices
K0e  = zeros(ndsPerElem*dofN,ndsPerElem*dofN); 
% B = zeros(dofN*dofN,ndsPerElem*dofN);

for i = 1:numel(int_point)
     
    xi = int_point(i);
   
    for j = 1:numel(int_point)
        
        niu = int_point(j);
        
        % Shape function derivatives (Natural Coords):
        R = 1/4.*[(-1+niu) (1-niu) (1+niu) (-1-niu); 
                   (-1+xi) (-1-xi)  (1+xi)   (1-xi)];
        
        % Jacobian Matrix
        J =  R * elCoords;
        
        % Shape function derivatives (Global Coords):
        dN = inv(J)*R;

        A1 = [2 0 0 0;
              0 2 0 0;
              0 0 1 1;
              0 0 1 1];
        
        % Strain Displacement Relationship
        B = B_Matrix_2D(dN);
     
        B0 = (1/2).*(A1 * B);
        
        K = (B0' * Cg * B0) .* det(J);
        
        K0e = K0e + int_weight(i)*int_weight(j).*K;

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
