function Pg = Assemble_Load_Vector(Q,input,dofPerNode,~,elementType)

% This script assembles the global load vector by generating the load
% vector for each given element. Stores in a sparse array

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% Q = displacement vector
% input = structure containing all input information for the problem

% OUTPUTS:
% Pg = sparse array storing global load vector at Q

%% ---------------------------------------------------------------------------
PROGRESS = 0;

% 
nNodes = size(input.ND,1);
totalDOFS = dofPerNode * nNodes;

Pg = zeros(totalDOFS,1);

% Number of elements
numE = size(input.EL,1);   

fibres = [];


% Set up progress bar for assembly
if PROGRESS == 1
    hwb    = waitbar(0,'Assembly process ...');
    set(hwb,'Name','Assembling Global Force Vector');
    indj = 1;
end

for ele=1:numE

    % Update progress bar
    if PROGRESS == 1
        if ele>(indj*50)
            waitbar(ele/numE,hwb,[num2str(floor(100*ele/numE)) ' % assembly ...']);
            indj=indj+1; 
        end
        if ele==numE
            waitbar(ele/numE,hwb,[' assembly finished ']);
        end
    end
    
    %%   2D QUADRILATERAL ELEMENT

    if elementType == 1
        
        % Equivalent global node number of local nodes
        node1 = find(input.ND(:,1)==input.EL(ele,3));
        node2 = find(input.ND(:,1)==input.EL(ele,4));
        node3 = find(input.ND(:,1)==input.EL(ele,5));
        node4 = find(input.ND(:,1)==input.EL(ele,6));
        
        % Local node coordinate info
        x1 = input.ND(node1,2); y1 = input.ND(node1,3);
        x2 = input.ND(node2,2); y2 = input.ND(node2,3);
        x3 = input.ND(node3,2); y3 = input.ND(node3,3);
        x4 = input.ND(node4,2); y4 = input.ND(node4,3);
        
        localCoords = [x1 y1; x2 y2; x3 y3; x4 y4];
        
        % Material orientation
        theta = input.EL(ele,end);
        
        % Local displacement vector
        id = [node1*dofPerNode-1; node1*dofPerNode-0;
              node2*dofPerNode-1; node2*dofPerNode-0;
              node3*dofPerNode-1; node3*dofPerNode-0;
              node4*dofPerNode-1; node4*dofPerNode-0];
        q = Q(id);
        
        %% Determine element load vector
        % Pe = Element_Load_Vector_2D(q,localCoords,theta,input.model);
        Pe = Element_Load_Vector_2D_TANGENT(q,localCoords,theta,input.model);
        
        %% Assemble element load vector to global load vector
        
        Pg(id) = Pg(id) + Pe;
       
    end

    if elementType == 2 || elementType == 4
        
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
        
        localCoords = [x1 y1 z1; 
                       x2 y2 z2; 
                       x3 y3 z3; 
                       x4 y4 z4; 
                       x5 y5 z5; 
                       x6 y6 z6; 
                       x7 y7 z7; 
                       x8 y8 z8];
        
        % 
        id = [node1*dofPerNode-2; node1*dofPerNode-1; node1*dofPerNode-0;
              node2*dofPerNode-2; node2*dofPerNode-1; node2*dofPerNode-0;
              node3*dofPerNode-2; node3*dofPerNode-1; node3*dofPerNode-0;
              node4*dofPerNode-2; node4*dofPerNode-1; node4*dofPerNode-0;
              node5*dofPerNode-2; node5*dofPerNode-1; node5*dofPerNode-0;
              node6*dofPerNode-2; node6*dofPerNode-1; node6*dofPerNode-0;
              node7*dofPerNode-2; node7*dofPerNode-1; node7*dofPerNode-0;
              node8*dofPerNode-2; node8*dofPerNode-1; node8*dofPerNode-0];
        
        % Displacements of local nodes
        q = Q(id);
        
        % Material orientation
        theta = input.EL(ele,end);
        
        %% Determine element load vector
%         Pe = Element_Load_Vector_3D(q,localCoords,theta,input.L0,input.M0);
        if isempty(input.FIBRES) == 0
            fibres.f = input.FIBRES.f(ele,:);
            fibres.s = input.FIBRES.s(ele,:);
            fibres.n = input.FIBRES.n(ele,:);
        end
        Pe = Element_Load_Vector_3D_TANGENT(q,localCoords,theta,input.model,fibres);
        
        %% Assemble element load vector to global load vector
        
        Pg(id) = Pg(id) + Pe;
       
    end
     
end

%% Store global stiffness matrix as a sparse matrix
% uses the vectors rowK, colK and valK to make a sparse matrix that is square and of dimension #nodes * DoF_per_node
% Kgl(idx) = valK(rowK(idx), colK(idx)) is how this is generated

% Pg = sparse(rowK,colK,valK,dofN*size(input.ND,1),1);

% AA = full(Kgl)      % use to view check completed global stiffness matrix
if PROGRESS == 1
    close(hwb)
end

end