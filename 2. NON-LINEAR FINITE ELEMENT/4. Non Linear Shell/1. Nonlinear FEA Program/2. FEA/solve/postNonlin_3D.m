function [F_,C_,E0_,En_,S_,tau_,sig_,P_,q_] = postNonlin_3D(input,solvedQ,solvedP,eleAnalyse)

% This script runs the post-processing for a 3D problem in the nonlinear
% finite element program. It returns arrays storing the deformation
% measures and stress tensors at each node + centre of each element

% Author:   Daniel O'Shea
% Created:  13 April 2018

% INPUTS:
% input = contains all input information
% Qsol = solved displacement field at each load step
% Psol = load vector at each load step

% OUTPUTS:
% F_ = deformation gradient tensor
% C_ = right Cauchy-Green deformation tensor
% S_ = Second Piola Kirchhoff stress tensor
% tau_ = Kirchhoff stress tensor


%% -----------------------------------------------------------------------

numElem = numel(input.EL(eleAnalyse,1));

dofPerNode = 3;
nComponents = dofPerNode * dofPerNode;
nodesPerElem = 8;
loadSteps = size(solvedP,2);

% [node1 node2 node3 node4 centre] 
pointEta = [ -1  -1   1   1  -1 -1  1  1  0]; 
pointXi  = [ -1   1   1  -1  -1  1  1 -1  0]; 
pointMu  = [ -1  -1  -1  -1   1  1  1  1  0];

% Initialise desired kinematic and stress tensors
F_ = zeros(numElem,(nodesPerElem+1)*nComponents,loadSteps);
C_ = zeros(numElem,(nodesPerElem+1)*nComponents,loadSteps);
S_ = zeros(numElem,(nodesPerElem+1)*nComponents,loadSteps);
P_ = zeros(numElem,(nodesPerElem+1)*nComponents,loadSteps);
tau_ = zeros(numElem,(nodesPerElem+1)*nComponents,loadSteps);
sig_ = zeros(numElem,(nodesPerElem+1)*nComponents,loadSteps);
E0_ = 1;    % need to complete
En_ = 1;    % need to complete
q_ = zeros(numElem,nodesPerElem*dofPerNode,loadSteps);

% Identity tensor
I = eye(3);

fibres = [];

% for each load step
for load = 1 : loadSteps

    Q = solvedQ(:,load);
    
    % for each element
    for ele = 1:numElem

        % Current element material orientation
        theta = input.EL(ele,end);
        phi = deg2rad(theta);

        if ~isempty(input.FIBRES)
            fibres.f = input.FIBRES.f(ele,:);
            fibres.s = input.FIBRES.s(ele,:);
            fibres.n = input.FIBRES.n(ele,:);
        end

        % Global Element Constitutive Law
        % [~,~,T,~,~]=Transform2(phi);
        % Lambda = inv(T) * Lambda * T ;
        % Mu     = inv(T) *   Mu   * T ;
        
         % Equivalent global node number of local nodes
        node1 = find(input.ND(:,1)==input.EL(ele,3));
        node2 = find(input.ND(:,1)==input.EL(ele,4));
        node3 = find(input.ND(:,1)==input.EL(ele,5));
        node4 = find(input.ND(:,1)==input.EL(ele,6));
        node5 = find(input.ND(:,1)==input.EL(ele,7));
        node6 = find(input.ND(:,1)==input.EL(ele,8));
        node7 = find(input.ND(:,1)==input.EL(ele,9));
        node8 = find(input.ND(:,1)==input.EL(ele,10));
        
        % Local node coordinate info
        x1 = input.ND(node1,2); y1 = input.ND(node1,3); z1 = input.ND(node1,4);
        x2 = input.ND(node2,2); y2 = input.ND(node2,3); z2 = input.ND(node2,4);
        x3 = input.ND(node3,2); y3 = input.ND(node3,3); z3 = input.ND(node3,4);
        x4 = input.ND(node4,2); y4 = input.ND(node4,3); z4 = input.ND(node4,4);
        x5 = input.ND(node5,2); y5 = input.ND(node5,3); z5 = input.ND(node5,4);
        x6 = input.ND(node6,2); y6 = input.ND(node6,3); z6 = input.ND(node6,4);
        x7 = input.ND(node7,2); y7 = input.ND(node7,3); z7 = input.ND(node7,4);
        x8 = input.ND(node8,2); y8 = input.ND(node8,3); z8 = input.ND(node8,4);
        
        elCoords = [x1 y1 z1; 
                    x2 y2 z2; 
                    x3 y3 z3; 
                    x4 y4 z4; 
                    x5 y5 z5; 
                    x6 y6 z6; 
                    x7 y7 z7; 
                    x8 y8 z8];
        
        id = [node1*dofPerNode-2; node1*dofPerNode-1; node1*dofPerNode-0;
              node2*dofPerNode-2; node2*dofPerNode-1; node2*dofPerNode-0;
              node3*dofPerNode-2; node3*dofPerNode-1; node3*dofPerNode-0;
              node4*dofPerNode-2; node4*dofPerNode-1; node4*dofPerNode-0;
              node5*dofPerNode-2; node5*dofPerNode-1; node5*dofPerNode-0;
              node6*dofPerNode-2; node6*dofPerNode-1; node6*dofPerNode-0;
              node7*dofPerNode-2; node7*dofPerNode-1; node7*dofPerNode-0;
              node8*dofPerNode-2; node8*dofPerNode-1; node8*dofPerNode-0];
        q = Q(id);
        q_(ele,:,load) = q';

        F_store = zeros(1,(nodesPerElem+1)*nComponents);
        C_store = zeros(1,(nodesPerElem+1)*nComponents);
        S_store = zeros(1,(nodesPerElem+1)*nComponents);
        P_store = zeros(1,(nodesPerElem+1)*nComponents);
        tau_store = zeros(1,dofPerNode*nComponents);
        
        for j = 1:numel(pointXi)
            
            if load == 63
                hello = 0;
            end

            % Current node location in natural coordinates
            xi = pointXi(j);
            eta = pointEta(j);
            mu  = pointMu(j);
            
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
            % A1m = T2M(A1);

            % Get PK2 stress tensor, and tangent stiffness tensor           
            [Smat,~] = Get_S_CC(F, input.model, fibres);
            Svec = M2V(Smat);

            % Deformation Gradient vector
            Fvec = B * q + [1;1;1;0;0;0;0;0;0];
            Fmat = V2M(Fvec);
            Cvec = M2V(F.'*F);

            % Kirchoff Stresses:
            tauMat = Fmat*Smat*Fmat.';
            tauVec = M2V(tauMat);

            % Cauchy Stresses:
            sigVec = tauVec ./ det(F);
            
            if ~isreal(sigVec)
                hello = 64;
            end

            % First Piola Kirchoff Stresses:
            Pmat = Fmat * transpose(Smat);
            Pvec = M2V(Pmat);

            % Store results at all nodes
            F_store(1,j*nComponents-(nComponents-1):j*nComponents) = Fvec';
            C_store(1,j*nComponents-(nComponents-1):j*nComponents) = Cvec';
            S_store(1,j*nComponents-(nComponents-1):j*nComponents) = Svec';
            P_store(1,j*nComponents-(nComponents-1):j*nComponents) = Pvec';            
            tau_store(1,j*nComponents-(nComponents-1):j*nComponents) = tauVec';
            sig_store(1,j*nComponents-(nComponents-1):j*nComponents) = sigVec';
            
        end
        
        F_(ele,:,load) = F_store;
        C_(ele,:,load) = C_store;
        S_(ele,:,load) = S_store;
        P_(ele,:,load) = P_store;
        tau_(ele,:,load) = tau_store;
        sig_(ele,:,load) = sig_store;
        
        
        
    end
    
        
end