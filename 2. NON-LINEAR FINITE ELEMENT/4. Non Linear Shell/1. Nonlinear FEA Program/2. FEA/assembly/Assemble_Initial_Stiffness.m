function K0 = Assemble_Initial_Stiffness(input,type,ndsE,dofN)

% This script assembles the global stiffness matrix for the infinitesimal
% case

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% input = strucure containing all input information for the problem
% type = element type
% ndsPerElem = nodes per element
% dofN = degrees of freedom per node

% OUTPUTS:
% K0 = sparse array storing initial global stiffness matrix

%% ---------------------------------------------------------------------------
numE = size(input.EL,1);      % number of elements
numN = size(input.ND,1);      % number of nodes
dofT = dofN * numN;

% Total Degrees of Freedom per element
dofE = ndsE*dofN;

% matrix size for a single element    
Lsq   = dofE * dofE;
    
% total DoF for the entire system
lngth = Lsq * numE;
    
% initalise global stiffness matrix identifiers
rowK  = zeros(lngth,1);
colK  = zeros(lngth,1);
valK  = zeros(lngth,1);
cIdx  = 1;

fibres = [];

for i=1:numE
    
    %%   2D QUADRILATERAL ELEMENT

    if type == 1 || type == 3
        
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
        
        localCoords = [x1 y1; x2 y2; x3 y3; x4 y4];
        
        % Material orientation
        theta = input.EL(i,end);

        % Determine element stiffness matrix
        K0e = Element_Initial_Stiffness_2D(localCoords,theta,input.model);
                
        % Assemble element stiffness to global stiffness
        t1 = node1*dofN; t2 = node2*dofN; t3 = node3*dofN; t4 = node4*dofN;
        bg = (t1-(dofN-1)):t1;  md1= (t2-(dofN-1)):t2;
        md2= (t3-(dofN-1)):t3;  en = (t4-(dofN-1)):t4;

        qL  = [bg md1 md2 en];

        rowLoc = reshape(repmat(qL,dofE,1),Lsq,1);
        colLoc = repmat(qL',dofE,1);

        linAr = cIdx:cIdx+Lsq-1;
        rowK(linAr) = rowK(linAr) + rowLoc;
        colK(linAr) = colK(linAr) + colLoc;
        valK(linAr) = valK(linAr) + K0e(:);

        cIdx = cIdx + Lsq;
       
    end

    %%   3D BRICK ELEMENT
    
    if type == 2 || type == 4
        
        % Equivalent global node number of local nodes
        node1 = find(input.ND(:,1)==input.EL(i,3));
        node2 = find(input.ND(:,1)==input.EL(i,4));
        node3 = find(input.ND(:,1)==input.EL(i,5));
        node4 = find(input.ND(:,1)==input.EL(i,6));
        node5 = find(input.ND(:,1)==input.EL(i,7));
        node6 = find(input.ND(:,1)==input.EL(i,8));
        node7 = find(input.ND(:,1)==input.EL(i,9));
        node8 = find(input.ND(:,1)==input.EL(i,10));
        
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
        
        % Material orientation
        theta = input.EL(i,end);

        % Determine element stiffness matrix
        if isempty(input.FIBRES) == 0
            fibres.f = input.FIBRES.f(i,:);
            fibres.s = input.FIBRES.s(i,:);
            fibres.n = input.FIBRES.n(i,:);
        end
        K0e = Element_Initial_Stiffness_3D(localCoords,theta,input.model, fibres);
                
        % Assemble element stiffness to global stiffness
        t1 = node1*dofN; t2 = node2*dofN; t3 = node3*dofN; t4 = node4*dofN;
        t5 = node5*dofN; t6 = node6*dofN; t7 = node7*dofN; t8 = node8*dofN;

        bg = [(t1-(dofN-1)):t1];
        md1= [(t2-(dofN-1)):t2];
        md2= [(t3-(dofN-1)):t3];
        md3= [(t4-(dofN-1)):t4];
        md4= [(t5-(dofN-1)):t5];
        md5= [(t6-(dofN-1)):t6];
        md6= [(t7-(dofN-1)):t7];
        en = [(t8-(dofN-1)):t8];
        
        qL     = [bg md1 md2 md3 md4 md5 md6 en];
        rowLoc = reshape(repmat(qL,dofE,1),Lsq,1); 
        colLoc = repmat(qL',dofE,1);  
        linAr  = cIdx:cIdx+Lsq-1;
        rowK(linAr) = rowK(linAr) + rowLoc;
        colK(linAr) = colK(linAr) + colLoc;
        valK(linAr) = valK(linAr) + K0e(:);
        
        cIdx = cIdx + Lsq; 
       
    end
    
end

%% Store global stiffness matrix as a sparse matrix
% uses the vectors rowK, colK and valK to make a sparse matrix that is square and of dimension #nodes * DoF_per_node
% Kgl(idx) = valK(rowK(idx), colK(idx)) is how this is generated

K0 = sparse(rowK,colK,valK,dofT,dofT);

% Close the progress bar

end