function [Kt, thetaElems, hElems] = Assemble_Tangent_Stiffness_MIXED(Q, input, type, ndPerElem, dofN, PROG, lambdaElems, thetaElems, hElems)

% This script assembles the tangent stiffness matrix for at a given
% displacement field, using the absolute tensor expression derived by Mario
% Attard

% Author:   Daniel O'Shea
% Created:  22 March 2022

% INPUTS:
% Q = current global displacement vector
% input = strucure containing all input information for the problem
% type = element type
% ndsPerElem = nodes per element
% dofN = degrees of freedom per node

% OUTPUTS:
% Ks = sparse array storing secant global stiffness matrix

% UPDATE LOG:


%% ---------------------------------------------------------------------------

% Total Degrees of Freedom per element
dofE = ndPerElem*dofN;

% matrix size for a single element    
Lsq   = dofE*dofE;
    
% total DoF for the entire system
lngth = (dofE)^2*size(input.EL,1);
    
% initalise global stiffness matrix identifiers
rowK  = zeros(lngth,1);
colK  = zeros(lngth,1);
valK  = zeros(lngth,1);
cIdx  = 1;

numElem = size(input.EL,1);      % number of elements

% Set up progress bar for assembly
if PROG == 1
    hwb    = waitbar(0,'Assembly process ...');
    set(hwb,'Name','Assembling Global Tangent Stiffness Matrix');
    indj = 1;
end

for i=1:numElem
    
    % Update progress bar
    if PROG == 1
        if i>(indj*50)
            waitbar(i/numElem,hwb,[num2str(floor(100*i/numElem)) ' % assembly ...']);
            indj=indj+1; 
        end
        if i==numElem
            waitbar(i/numElem,hwb,[' assembly finished ']);
    %         disp(['     K assembly: ' num2str(i) ' elements - complete']);  
        end
    end

    %%   2D QUADRILATERAL ELEMENT

    if type == 1
        
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
        id = [node1*dofN-1;
              node1*dofN-0;
              node2*dofN-1;
              node2*dofN-0;
              node3*dofN-1;
              node3*dofN-0;
              node4*dofN-1;
              node4*dofN-0];
        q = Q(id);
        
        % Material orientation
        theta = input.EL(i,end);

        lambdaE = lambdaElems(i);

        % Determine element stiffness matrix
        mixed = 1;
        if mixed == 1
            % we want to pass in 'current coords' as localCoords here...
            [Ke, thetaE] = Element_Tangent_Stiffness_2D_MIXED(q,localCoords,theta,input.model,lambdaE);
            thetaElems(i) = thetaE;
            hElems(i) = (thetaE^2 - 1)/thetaE;  % can adjust this function
        else
            Ke = Element_Tangent_Stiffness_2D(q,localCoords,theta,input.model);
        end
                
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
        valK(linAr) = valK(linAr) + Ke(:);

        cIdx = cIdx + Lsq;
       
    end

    %%  3D BRICK

    if type == 2
        
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
        
        id = [node1*dofN-2; node1*dofN-1; node1*dofN-0;
              node2*dofN-2; node2*dofN-1; node2*dofN-0;
              node3*dofN-2; node3*dofN-1; node3*dofN-0;
              node4*dofN-2; node4*dofN-1; node4*dofN-0;
              node5*dofN-2; node5*dofN-1; node5*dofN-0;
              node6*dofN-2; node6*dofN-1; node6*dofN-0;
              node7*dofN-2; node7*dofN-1; node7*dofN-0;
              node8*dofN-2; node8*dofN-1; node8*dofN-0];
        q = Q(id);
        
        % Material orientation
        theta = input.EL(i,end);

        lambdaE = lambdaElems(i);

        % Determine element stiffness matrix
        mixed = 1;
        if mixed == 1
            % we want to pass in 'current coords' as localCoords here...
            [Ke, thetaE] = Element_Tangent_Stiffness_3D_MIXED(q,localCoords,theta,input.model,lambdaE);
            thetaElems(i) = thetaE;
            hElems(i) = (thetaE^2 - 1)/thetaE; % can adjust this function
        else
            Ke = Element_Tangent_Stiffness_3D(q,localCoords,theta,input.model);
        end

        % Assemble element stiffness to global stiffness
        t1 = node1*dofN; t2 = node2*dofN; t3 = node3*dofN; t4 = node4*dofN;
        t5 = node5*dofN; t6 = node6*dofN; t7 = node7*dofN; t8 = node8*dofN;

        bg = (t1-(dofN-1)):t1;
        md1= (t2-(dofN-1)):t2;
        md2= (t3-(dofN-1)):t3;
        md3= (t4-(dofN-1)):t4;
        md4= (t5-(dofN-1)):t5;
        md5= (t6-(dofN-1)):t6;
        md6= (t7-(dofN-1)):t7;
        en = (t8-(dofN-1)):t8;
        
        qL     = [bg md1 md2 md3 md4 md5 md6 en];
        rowLoc = reshape(repmat(qL,dofE,1),Lsq,1); 
        colLoc = repmat(qL',dofE,1);  
        linAr  = cIdx:cIdx+Lsq-1;
        rowK(linAr) = rowK(linAr) + rowLoc;
        colK(linAr) = colK(linAr) + colLoc;
        valK(linAr) = valK(linAr) + Ke(:);
        
        cIdx = cIdx + Lsq; 
       
    end
    
end

%% Store global stiffness matrix as a sparse matrix
% uses the vectors rowK, colK and valK to make a sparse matrix that is square and of dimension #nodes * DoF_per_node
% Kgl(idx) = valK(rowK(idx), colK(idx)) is how this is generated

Kt = sparse(rowK,colK,valK,dofN*size(input.ND,1),dofN*size(input.ND,1));

% AA = full(Kgl)      % use to view check completed global stiffness matrix

% Close the progress bar
if PROG == 1
    close(hwb);
end

end