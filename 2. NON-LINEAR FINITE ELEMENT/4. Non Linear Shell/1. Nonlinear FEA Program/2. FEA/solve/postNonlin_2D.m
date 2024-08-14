function [F_,C_,E0_,En_,S_,tau_,P_,q_] = postNonlin_2D(input,Qsol,Psol,ele_analyse)

% This script runs the post-processing for a 2D problem in the nonlinear
% finite element program. It returns arrays storing the deformation
% measures and stress tensors at each node + centre of each element

% Author:   Daniel O'Shea
% Created:  23 March 2018

% INPUTS:
% input = contains all input information
% Qsol = solved displacement field at each load step
% Psol = load vector at each load step

% OUTPUTS:
% F_ = deformation gradient tensor
% C_ = right Cauchy-Green deformation tensor
% E0_ = logarithmic strain tensor
% En_ = polynomial strain tensor
% S_ = Second Piola Kirchhoff stress tensor
% tau_ = Kirchhoff stress tensor
% q_ = local element displacement vector

% UPDATE LOG
% 26/04/2018 - Added strain tensors to analysis and output
% 26/04/2018 - Added local displacement vector as output

%% -----------------------------------------------------------------------

numElem = numel(input.EL(ele_analyse,1));
numNode = numel(input.ND(:,1));

dofN = 2;
numGrad = dofN * dofN;
ndsPerElem = 4;
N = size(Psol,2);

I = eye(2);

% [node1 node2 node3 node4 centre] 
point_eta = [ -1  -1   1   1  0]; 
point_xi  = [ -1   1   1  -1  0]; 

F_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
C_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
E0_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
En_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
S_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
tau_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
P_ = zeros(numElem,(ndsPerElem+1)*numGrad,N);
q_ = zeros(numElem,ndsPerElem*dofN,N);

% for each load step
for nn = 1 : N

    Q = Qsol(:,nn);
    
    % for each element
    for i = 1:numElem
        
        % Principle Element Constitutive Law
        % Lambda = input.L0;
        % Mu = input.M0;
        
        % Current element material orientation
        theta = input.EL(i,end);
        phi = deg2rad(theta);

        % Global Element Constitutive Law
        [T4,~,~,~,~]=Transform2(phi);
        % Lambda = inv(T4) * Lambda * T4 ;
        % Mu     = inv(T4) *   Mu   * T4 ;
        
        % Equivalent global node number of local nodes
        node1 = find(input.ND(:,1)==input.EL(i,3));
        node2 = find(input.ND(:,1)==input.EL(i,4));
        node3 = find(input.ND(:,1)==input.EL(i,5));
        node4 = find(input.ND(:,1)==input.EL(i,6));

        % Local node coordinate info
        x1 = input.ND(node1,2); y1 = input.ND(node1,3);
        x2 = input.ND(node2,2); y2 = input.ND(node2,3);
        x3 = input.ND(node3,2); y3 = input.ND(node3,3);
        x4 = input.ND(node4,2); y4 = input.ND(node4,3);

        elCoords = [x1 y1; x2 y2; x3 y3; x4 y4];
        id = [node1*dofN-1;
              node1*dofN-0;
              node2*dofN-1;
              node2*dofN-0;
              node3*dofN-1;
              node3*dofN-0;
              node4*dofN-1;
              node4*dofN-0];
        q = Q(id);
        q_(i,:,nn) = q';
        
        F_store = zeros(1,(ndsPerElem+1)*numGrad);
        C_store = zeros(1,(ndsPerElem+1)*numGrad);
        E0_store = zeros(1,(ndsPerElem+1)*numGrad);
        En_store = zeros(1,(ndsPerElem+1)*numGrad);
        S_store = zeros(1,(ndsPerElem+1)*numGrad);
        P_store = zeros(1,(ndsPerElem+1)*numGrad);
        tau_store = zeros(1,dofN*numGrad);
        for j = 1:numel(point_xi)
            xi = point_xi(j);
            eta = point_eta(j);

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
            % A1m = T2M(A1);

            % Get PK2 stress tensor, and tangent stiffness tensor           
            [Smat,~] = Get_S_CC(GradU, input.model);
            Svec = M2V(Smat);

            % Deformation Gradient vector
            Fvec = B * q + [1;1;0;0];
            Fmat = V2M(Fvec);
            Cvec = M2V(F.'*F);
            
            % Kirchoff Stresses:
            taumat = transpose(Fmat)*Smat*Fmat;
            tauvec = M2V(taumat);
            
            % First Piola Kirchoff Stresses:
            Pmat = Fmat * transpose(Smat);
            Pvec = M2V(Pmat);

            % Store results at all nodes
            F_store(1,j*numGrad-(numGrad-1):j*numGrad) = Fvec';
            C_store(1,j*numGrad-(numGrad-1):j*numGrad) = Cvec';
            % E0_store(1,j*numGrad-(numGrad-1):j*numGrad)= E0vec';
            % En_store(1,j*numGrad-(numGrad-1):j*numGrad)= Envec(:,1)';
            S_store(1,j*numGrad-(numGrad-1):j*numGrad) = Svec';
            P_store(1,j*numGrad-(numGrad-1):j*numGrad) = Pvec';
            tau_store(1,j*numGrad-(numGrad-1):j*numGrad) = tauvec';
        end
        
        F_(i,:,nn) = F_store;
        C_(i,:,nn) = C_store;
        % E0_(i,:,nn)= E0_store;
        % En_(i,:,nn)= En_store;
        S_(i,:,nn) = S_store;
        P_(i,:,nn) = P_store;
        tau_(i,:,nn) = tau_store;

        E0_ = 1;    % need to complete
        En_ = 1;    % need to complete

    end
    
        
end