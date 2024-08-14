%% plotStressContour

% This script produces contour plots of the stress field for a 2D nonlinear
% finite element

% Author:   Daniel O'Shea
% Created:  22 March 2018


%% ------------------------------------------------------------------------
numNode = size(input.ND,1);
numElem = size(input.EL,1);
dofN = 2;

point_eta = [ -1  -1   1   1]; 
point_xi  = [ -1   1   1  -1]; 

Stress_store = zeros(numElem,1+4*4);
Stress_store(:,1) = 1:1:numElem;
Lambda = input.L0;
Mu = input.M0;

for nn = 1 : N
    figure
    Q = Qsol(:,nn);
    
    for i = 1:numElem

        % Current element material orientation
        theta = input.EL(i,end);
        phi = deg2rad(theta);

        % Global Element Constitutive Law
        [T4,~,~,~,~]=Transform2(phi);
        Lambda = inv(T4) * Lambda * T4 ;
        Mu     = inv(T4) *   Mu   * T4 ;
        
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

        stress = zeros(4,4);
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

            % Displacement Gradient Vector
            GradU = B * q;

            % Assemble A1 matrix:
            A1 = [2 0 0 0;
                  0 2 0 0;
                  0 0 1 1;
                  0 0 1 1];

            % Assemble A2 matrix:
            A2 = [GradU(1)   0        0      GradU(4);
                      0    GradU(2) GradU(3)   0     ;
                      0    GradU(4) GradU(1)   0     ;
                  GradU(3)   0        0      GradU(2)];

            % Assemble A3 matrix:
            A3 = [2*GradU(1)     0          0      2*GradU(4);
                      0      2*GradU(2) 2*GradU(3)     0     ;
                  1*GradU(3) 1*GradU(4) 1*GradU(1) 1*GradU(2);
                  1*GradU(3) 1*GradU(4) 1*GradU(1) 1*GradU(2)];

            % Right Cauchy-Green deformation vector:
            Cvec = [1;1;0;0] + (A1 + A2)*GradU;

            % Compute stress vector:
            Svec = S_Vector_2D(Cvec,Lambda,Mu);
            Smat = V2M(Svec);
            
            % Deformation Gradient vector
            Fvec = B * q + [1;1;0;0];
            Fmat = V2M(Fvec);
            
            % Kirchoff Stresses:
            taumat = transpose(Fmat)*Smat*Fmat;
            tauvec = M2V(taumat);
            
            stress(:,j) = tauvec;
            
        end

        Stress_store(i,2:end,nn) = reshape(stress,1,4*4);
        str11 = reshape(stress(2,:),2,2)';
        X = [x1 x2; x4 x3];
        Y = [y1 y2; y4 y3];
        contourf(X,Y,str11)
        colorbar
        hold on
        caxis('manual')
    end
    
    
end