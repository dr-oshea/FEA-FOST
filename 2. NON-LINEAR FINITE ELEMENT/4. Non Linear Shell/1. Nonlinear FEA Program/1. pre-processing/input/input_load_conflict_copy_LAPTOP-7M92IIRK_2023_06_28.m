function load = input_load(coords,element,xdim,ydim,zdim,mdiv,type,Force,Stress,FCASE)

% This function generates an array of nodes on which nodal displacements
% are applied, specifies the magnitude of that load and in which degree of
% freedom it acts

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% coords = nodal coordinate data
% xdim = total magnitude of geometry along x-dimension
% ydim = total magnitude of geometry along y-dimension
% type = element type
% Force = Applied nodal force [N]
% FCASE = force case

% OUTPUTS:
% load = array storing constraint information for applicable nodes

%% ---------------------------------------------------------------------------
tol = 10^-6; % tolerance for finding nodes

%% For a CSTQ element
if type == 1
    
    xmin = min(coords(:,2));
    ymin = min(coords(:,3));
    
    spacing = (xdim-xmin)/mdiv;
    
    % Significant node locations:
    
    % Corners
    p_lcorner = find(abs(coords(:,2)-xmin)<=tol&abs(coords(:,3)-ydim)<=tol);
    p_rcorner = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ydim)<=tol);
    n_lcorner = find(abs(coords(:,2)-xmin)<=tol&abs(coords(:,3)-ymin)<=tol);
    n_rcorner = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ymin)<=tol);
    
    % Faces/Edges:
    n_yface = find(abs(coords(:,3)-ymin)<=tol);
    p_yface = find(abs(coords(:,3)-ydim)<=tol);
    n_xface = find(abs(coords(:,2)-xmin)<=tol);
    p_xface = find(abs(coords(:,2)-xdim)<=tol);
    
    % Centre of Mesh
    centnode = find(coords(:,2)==xdim/2&coords(:,3)==ydim/2);
    
    
    %% Loads
    F_edge = Stress * spacing * 1;  % 2D
    F_corn = F_edge/2;
    
    % Pinned bottom edge, tension top edge in y-direction
    if FCASE == 1    
        load = zeros(numel(p_yface),3);
        for i = 1:numel(p_yface)    % top nodes
            if p_yface(i) == p_lcorner || p_yface(i) == p_rcorner
                load(i,:)=[p_yface(i) 0 F_corn];
            else
                load(i,:)=[p_yface(i) 0 F_edge];
            end
        end
               
    end
    
    % Uniaxial Tension
    if FCASE == 2
        
        load = zeros(numel(p_xface)+numel(n_xface),3);
        for i = 1:numel(p_xface)    % right nodes
            if p_xface(i) == p_rcorner || p_xface(i) == n_rcorner
                load(i,:)=[p_xface(i) F_corn 0];
            else
                load(i,:)=[p_xface(i) F_edge 0];
            end
        end
        for j = 1:numel(n_xface)    % left nodes
            if n_xface(j) == p_lcorner || n_xface(j) == n_lcorner
                load(j+i,:)=[n_xface(j) -F_corn 0];
            else
                load(j+i,:)=[n_xface(j) -F_edge 0];
            end
        end       
    end
    
    % Ciarletta Shear
    if FCASE == 3
        
        load = zeros(numel(p_yface),3);
        for i = 1:numel(p_yface)    % top nodes
            if p_yface(i) == p_lcorner || p_yface(i) == p_rcorner
                load(i,:)=[p_yface(i) F_corn 0];
            else
                load(i,:)=[p_yface(i) F_edge 0];
            end
        end        
    end
    
    % Simple Shear
    if FCASE == 4
        F_edge = Pressure * spacing;
        F_corn = F_edge/2;
        
        load = zeros(numel(p_yface),3);
        for i = 1:numel(p_yface)    % top nodes
            if p_yface(i) == p_lcorner || p_yface(i) == p_rcorner
                load(i,:)=[p_yface(i) F_corn 0];
            else
                load(i,:)=[p_yface(i) F_edge 0];
            end
        end
        for j = 1:numel(n_yface)    % bottom nodes
            if n_yface(j) == n_lcorner || n_yface(j) == n_rcorner
                load(i+j,:)=[n_yface(j) -F_corn 0];
            else
                load(i+j,:)=[n_yface(j) -F_edge 0];
            end
        end
        for k = 1:numel(n_xface)    % left nodes
            if n_xface(k) == n_lcorner || n_xface(k) == p_lcorner
                load(i+j+k,:)=[n_xface(k) 0 -F_corn];
            else
                load(i+j+k,:)=[n_xface(k) 0 -F_edge];
            end
        end
        for l = 1:numel(p_xface)    % right nodes
            if p_xface(l) == n_rcorner || p_xface(l) == p_rcorner
                load(i+j+k+l,:)=[p_xface(l) 0 F_corn];
            else
                load(i+j+k+l,:)=[p_xface(l) 0 F_edge];
            end
        end
    end
    
    % Equibiaxial Tension
    if FCASE == 5
        
        load(1,:) = [p_lcorner -F_corn  F_corn];
        load(2,:) = [p_rcorner  F_corn  F_corn];
        load(3,:) = [n_lcorner -F_corn -F_corn];
        load(4,:) = [n_rcorner  F_corn -F_corn];
        
        a = 1; b = size(load,1);
        for i = 1:numel(n_yface)
            if n_yface(i) ~= n_lcorner && n_yface(i) ~= n_rcorner
                load(b+a,:) = [n_yface(i) 0 -F_edge];
                a = a + 1;
            end
        end
        a = 1; b = size(load,1);
        for i = 1:numel(p_yface)
            if p_yface(i) ~= p_lcorner && p_yface(i) ~= p_rcorner
                load(b+a,:) = [p_yface(i) 0 F_edge];
                a = a + 1;
            end
        end
        a = 1; b = size(load,1);
        for i = 1:numel(n_xface)
            if n_xface(i) ~= n_lcorner && n_xface(i) ~= p_lcorner
                load(b+a,:) = [n_xface(i) -F_edge 0];
                a = a + 1;
            end
        end
        a = 1; b = size(load,1);
        for i = 1:numel(p_xface)
            if p_xface(i) ~= n_rcorner && p_xface(i) ~= p_rcorner
                load(b+a,:) = [p_xface(i) F_edge 0];
                a = a + 1;
            end
        end    
    end
    
    % Uniaxial Tension 1 element
    if FCASE == 6
        load = [2 F_corn 0;
                4 F_corn 0];
    end
    
end

%% For a 3D Brick element
if type == 2
    
    xmin = min(coords(:,2));
    ymin = min(coords(:,3));
    zmin = min(coords(:,4));
    
    spacing = (xdim-xmin)/mdiv;
    
    % Significant node locations:    
    corn_xnypzn = find(abs(coords(:,2)-0)<=tol&abs(coords(:,3)-ydim)<=tol&abs(coords(:,4)-0)<=tol);
    corn_xpypzn = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ydim)<=tol&abs(coords(:,4)-0)<=tol);
    corn_xnypzp = find(abs(coords(:,2)-0)<=tol&abs(coords(:,3)-ydim)<=tol&abs(coords(:,4)-zdim)<=tol);
    corn_xpypzp = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ydim)<=tol&abs(coords(:,4)-zdim)<=tol);
    corn_xnynzn = find(abs(coords(:,2)-0)<=tol&abs(coords(:,3)-0)<=tol&abs(coords(:,4)-0)<=tol);
    corn_xpynzn = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-0)<=tol&abs(coords(:,4)-0)<=tol);
    corn_xnynzp = find(abs(coords(:,2)-0)<=tol&abs(coords(:,3)-0)<=tol&abs(coords(:,4)-zdim)<=tol);
    corn_xpynzp = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-0)<=tol&abs(coords(:,4)-zdim)<=tol);

    edge_xnyn = find(abs(coords(:,2)-0)<=tol&abs(coords(:,3)-0)<=tol);
    edge_xnyp = find(abs(coords(:,2)-0)<=tol&abs(coords(:,3)-ydim)<=tol);
    edge_xpyn = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-0)<=tol);
    edge_xpyp = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ydim)<=tol);

    edge_ynzn = find(abs(coords(:,3)-0)<=tol&abs(coords(:,4)-0)<=tol);
    edge_ynzp = find(abs(coords(:,3)-0)<=tol&abs(coords(:,4)-zdim)<=tol);
    edge_ypzn = find(abs(coords(:,3)-ydim)<=tol&abs(coords(:,4)-0)<=tol);
    edge_ypzp = find(abs(coords(:,3)-ydim)<=tol&abs(coords(:,4)-zdim)<=tol);
    
    edge_znxn = find(abs(coords(:,4)-0)<=tol&abs(coords(:,2)-0)<=tol);
    edge_znxp = find(abs(coords(:,4)-0)<=tol&abs(coords(:,2)-xdim)<=tol);
    edge_zpxn = find(abs(coords(:,4)-zdim)<=tol&abs(coords(:,2)-0)<=tol);
    edge_zpxp = find(abs(coords(:,4)-zdim)<=tol&abs(coords(:,2)-xdim)<=tol);
    
    face_xn = find(abs(coords(:,2)-0)<=tol);
    face_xp = find(abs(coords(:,2)-xdim)<=tol);
    face_yn = find(abs(coords(:,3)-0)<=tol);
    face_yp = find(abs(coords(:,3)-ydim)<=tol);
    face_zn = find(abs(coords(:,4)-0)<=tol);
    face_zp = find(abs(coords(:,4)-zdim)<=tol);
    face_xc = find(abs(coords(:,2)-xdim/2)<=tol);
    face_yc = find(abs(coords(:,3)-ydim/2)<=tol);
    face_zc = find(abs(coords(:,4)-zdim/2)<=tol);

    centnode = find(coords(:,2)==xdim/2&coords(:,3)==ydim/2&coords(:,4)==zdim/2);
    
    %% Loads
    
    % Uniaxial Tension, X-Direction
    if FCASE == 2
        
        F = Stress * xdim * zdim;
        fx = F/(mdiv*mdiv*(xdim/zdim));   % Total force on one element
        F_face = fx;
        F_edge = F_face/2;
        F_corn = F_edge/2;
    

        load = zeros(numel(face_xn)+numel(face_xp),4);
        
        for i = 1:numel(face_xp)    % positive x-nodes
            node = face_xp(i);
            if ismember(node,[corn_xpynzn corn_xpynzp corn_xpypzn corn_xpypzp]) == 1
                load(i,:)=[face_xp(i) F_corn 0  0];
            elseif ismember(node,[edge_xpyn' edge_xpyp' edge_znxp' edge_zpxp']) == 1
                load(i,:)=[face_xp(i) F_edge 0  0];
            else
                load(i,:)=[face_xp(i) F_face 0  0];
            end
        end
        
        for j = 1:numel(face_xn)    % negative x-nodes
            node = face_xn(j);
            if ismember(node,[corn_xnynzn corn_xnynzp corn_xnypzn corn_xnypzp]) == 1
                load(i+j,:)=[face_xn(j) -F_corn 0  0];
            elseif ismember(node,[edge_xnyn' edge_xnyp' edge_znxn' edge_zpxn']) == 1
                load(i+j,:)=[face_xn(j) -F_edge 0  0];
            else
                load(i+j,:)=[face_xn(j) -F_face 0  0];
            end
        end
        
    end
        
    % Ciarletta Shear
    if FCASE == 3
        ndsPerElem = 8;
        
        F = Stress * xdim * zdim;       % Total force on face
        fx = F/(mdiv*mdiv*(xdim/zdim));   % Total force on one element
        fy = 0; fz = 0;
        
        topnodes = face_yp;
        load = zeros(numel(topnodes),3+1);
        load(:,1) = topnodes';

        % find elements attached to given node
        b=1;
        for i=1:numel(topnodes)
            node = topnodes(i);
            
            for j = 1:ndsPerElem
                el = find(element(:,j+2)==node);
                if isempty(el)==0; e(b)=el; b=b+1; end
            end
        end
        e=unique(e);
        
        % Apply forces to nodes in that element
        for j = 1:numel(e)
            el = e(j);
            index = find(ismember(element(el,3:3+ndsPerElem-1),topnodes));
            for k = 1:numel(index)
                % Row to store load info in array
                aa = find(load(:,1)==element(el,2+index(k)));
                
                if type == 2; factor = 1/4; end
%                 if EL_TYPE == 6 && ismember(index(k),1:4)==1; factor = 0; end
%                 if EL_TYPE == 6 && ismember(index(k),5:8)==1; factor = 0; end
                
                load(aa,2) = load(aa,2) + factor*fx;
                load(aa,3) = load(aa,3) + factor*fy;
                load(aa,4) = load(aa,4) + factor*fz;
            end
        end
    end
    
    % Equibiaxial Tension
    if FCASE == 5
        
        F = Stress * xdim * zdim;
        fx = F/(mdiv*mdiv*(xdim/zdim));   % Total force on one element
        F_face = fx;
        F_edge = F_face/2;
        F_corn = F_edge/2;
        
%         load(1,:) = [corn_xnypzn -F_corn  F_corn 0];
%         load(2,:) = [corn_xpypzn  F_corn  F_corn 0];
%         load(3,:) = [corn_xnypzp -F_corn  F_corn 0];
%         load(4,:) = [corn_xpypzp  F_corn  F_corn 0];
%         load(5,:) = [corn_xnynzn -F_corn -F_corn 0];
%         load(6,:) = [corn_xpynzn  F_corn -F_corn 0];
%         load(7,:) = [corn_xnynzp -F_corn -F_corn 0];
%         load(8,:) = [corn_xpynzp  F_corn -F_corn 0];
        
        % load = zeros(numel(face_xn)+numel(face_xp)+numel(face_yn)+numel(face_yp),4);
        load = zeros(numel(face_xp)+numel(face_yp),4);

        for i = 1:numel(face_xp)    % positive x-nodes
            node = face_xp(i);
            if ismember(node,[corn_xpynzn corn_xpynzp corn_xpypzn corn_xpypzp]) == 1
                load(i,:)=[face_xp(i) F_corn 0  0];
            elseif ismember(node,[edge_xpyn' edge_xpyp' edge_znxp' edge_zpxp']) == 1
                load(i,:)=[face_xp(i) F_edge 0  0];
            else
                load(i,:)=[face_xp(i) F_face 0  0];
            end
        end
        
        % for j = 1:numel(face_xn)    % negative x-nodes
        %     node = face_xn(j);
        %     if ismember(node,[corn_xnynzn corn_xnynzp corn_xnypzn corn_xnypzp]) == 1
        %         load(i+j,:)=[face_xn(j) -F_corn 0  0];
        %     elseif ismember(node,[edge_xnyn' edge_xnyp' edge_znxn' edge_zpxn']) == 1
        %         load(i+j,:)=[face_xn(j) -F_edge 0  0];
        %     else
        %         load(i+j,:)=[face_xn(j) -F_face 0  0];
        %     end
        % end
        % a = i+j;
        
        a = i;
        for i = 1:numel(face_yp)    % positive y-nodes
            node = face_yp(i);
            if ismember(node,[corn_xnypzn corn_xpypzn corn_xnypzp corn_xpypzp]) == 1
                load(a+i,:)=[face_yp(i) 0 F_corn 0];
            elseif ismember(node,[edge_xnyp' edge_xpyp' edge_ypzn' edge_ypzp']) == 1
                load(a+i,:)=[face_yp(i) 0 F_edge 0];
            else
                load(a+i,:)=[face_yp(i) 0 F_face 0];
            end
        end
        
        % for j = 1:numel(face_yn)    % negative y-nodes
        %     node = face_yn(j);
        %     if ismember(node,[corn_xnynzn corn_xpynzn corn_xnynzp corn_xpynzp]) == 1
        %         load(a+i+j,:)=[face_yn(j) 0 -F_corn 0 ];
        %     elseif ismember(node,[edge_xnyn' edge_xpyn' edge_ynzn' edge_ynzp']) == 1
        %         load(a+i+j,:)=[face_yn(j) 0 -F_edge 0 ];
        %     else
        %         load(a+i+j,:)=[face_yn(j) 0 -F_face 0 ];
        %     end
        % end
        
    end
    
    if FCASE == 6
        F_corn = 1;  % NEED TO FIX
        load = [2 F_corn 0 0;
                4 F_corn 0 0;
                6 F_corn 0 0;
                8 F_corn 0 0];
            
    end

    if FCASE == 91
        ndsPerElem = 8;
        
        F = Stress * xdim * zdim;       % Total force on face
        fx = F/(mdiv*mdiv*(xdim/zdim));   % Total force on one element
        fy = 0; fz = 0;
        
        topnodes = face_yp;
        load = zeros(numel(topnodes),3+1);
        load(:,1) = topnodes';

        % find elements attached to given node
        b=1;
        for i=1:numel(topnodes)
            node = topnodes(i);
            
            for j = 1:ndsPerElem
                el = find(element(:,j+2)==node);
                if isempty(el)==0; e(b)=el; b=b+1; end
            end
        end
        e=unique(e);
        
        % Apply forces to nodes in that element
        for j = 1:numel(e)
            el = e(j);
            index = find(ismember(element(el,3:3+ndsPerElem-1),topnodes));
            for k = 1:numel(index)
                % Row to store load info in array
                aa = find(load(:,1)==element(el,2+index(k)));
                
                if type == 2; factor = 1/4; end
%                 if EL_TYPE == 6 && ismember(index(k),1:4)==1; factor = 0; end
%                 if EL_TYPE == 6 && ismember(index(k),5:8)==1; factor = 0; end
                
                load(aa,2) = load(aa,2) + factor*fx;
                load(aa,3) = load(aa,3) + factor*fy;
                load(aa,4) = load(aa,4) + factor*fz;
            end
        end
        % Load case for FS mode
        % What we need is similar to FCASE = 3 above (though recreate here)
        % There is probably a ore efficient way to do it.. I made the code
        % in lines 354-363 more recently which achieves a similar result
        % though more efficiently

        % Note that we are defining a 'stress' being applied at nodes. To
        % convert to forces acting on nodes, we multiply stress by the area
        % that each node is 'servicing'. This means the force applied to a
        % node on the edge of a face, or at a vertex of cube will be less
        % that one applied to a node at the middle of a face (in fact, edge
        % = face / 2, corner = face / 4
        
        %load = [];
        % return 'load' array containing nodes with loads attached
    end

    if FCASE == 92
        % Load case for FN mode

    end
    
    if FCASE == 93
        % Load case for SF mode

    end
    
    if FCASE == 94
        % Load case for SN mode

    end
    
    if FCASE == 95
        % Load case for NF mode

    end
    
    if FCASE == 96
        % Load case for NS mode

    end

end

end