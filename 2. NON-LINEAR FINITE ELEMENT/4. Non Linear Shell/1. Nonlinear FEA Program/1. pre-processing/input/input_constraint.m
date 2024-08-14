function cons = input_constraint(coords,xdim,ydim,zdim,type,CCASE,mesh_div,aspect,thickness_div,height_div)

% This function generates an array of nodes on which displacement is
% constrained, and specifies which degree of freedom at that node

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% coords = nodal coordinate data
% xdim = total magnitude of geometry along x-dimension
% ydim = total magnitude of geometry along y-dimension
% type = element type
% CCASE = constrain case

% OUTPUTS:
% cons = array storing constraint information for applicable nodes

%% ---------------------------------------------------------------------------
tol = 10^-6; % tolerance for finding nodes

%% For a CSTQ element (2D)
if type == 1 || type == 3
    
    xmin = min(coords(:,2));
    ymin = min(coords(:,3));
    
    % Significant node locations:
    
    % Corners
    p_lcorner = find(abs(coords(:,2)-xmin)<=tol&abs(coords(:,3)-ydim)<=tol);
    p_rcorner = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ydim)<=tol);
    n_lcorner = find(abs(coords(:,2)-xmin)<=tol&abs(coords(:,3)-ymin)<=tol);
    n_rcorner = find(abs(coords(:,2)-xdim)<=tol&abs(coords(:,3)-ymin)<=tol);
    
    % Faces/Edges:
    n_yface = find(abs(coords(:,3)-ymin)<=tol);
    c_yface = find(abs(coords(:,3)-ydim/2)<=tol);
    p_yface = find(abs(coords(:,3)-ydim)<=tol);
    n_xface = find(abs(coords(:,2)-xmin)<=tol);
    c_xface = find(abs(coords(:,2)-xdim/2)<=tol);
    p_xface = find(abs(coords(:,2)-xdim)<=tol);
    
    % Centre of Mesh
    centnode = find(coords(:,2)==xdim/2&coords(:,3)==ydim/2);
    
    
    %% Constraints
    
    % 2DOF per node, pinned bottom edge
    if CCASE == 1
        cons = zeros(numel(n_yface),3);
        for i = 1:numel(n_yface)    % bottom nodes
            cons(i,:)=[n_yface(i) 0 0];
        end
    end
    
    % 2DOF per node, centre x-row nodes constrained
    if CCASE == 2
        cons = zeros(numel(c_xface),3);
        for i = 1:numel(c_xface)   % centre nodes
            if c_xface(i) == centnode
                cons(i,:) = [c_xface(i) 0 0];
            else
                cons(i,:) = [c_xface(i) 0 1];
            end
        end
    end
    
    % 2DOF per node, Ciarletta shear (top edge roller, botom edge pin)
    if CCASE == 3
        cons = zeros(numel(n_yface),3);
        for i = 1:numel(n_yface)    % bottom nodes
            cons(i,:)=[n_yface(i) 0 0];
        end
        for j = 1:numel(p_yface)    % top nodes
            cons(i+j,:)=[p_yface(j) 1 1];
        end
    end
    
    % 2DOF per node, Simple shear
    if CCASE == 4
        % Need to complete
    end
    
    % 2DOF per node, Equibiaxial Tension
    if CCASE == 5
        cons = zeros(numel(c_xface),3);
        for i = 1:numel(c_xface)   % centre nodes
            if c_xface(i) == centnode
                cons(i,:) = [c_xface(i) 0 0];
            else
                cons(i,:) = [c_xface(i) 0 1];
            end
        end
%         for j = 1:numel(c_yface)
%             if c_yface(j) ~= centnode
%                 cons(i+j,:) = [c_yface(j) 1 0];
%             end
%         end
    end
    
    % 2DOF per node, Uniaixal tension for 1 element
    if CCASE == 6
        cons = [1 0 0;
                2 1 0;
                3 0 1];
    end
    
    % 2DOF per node, Pure Shear Biaxial Test
    if CCASE == 7
        for i = 1:numel(c_xface)
            if c_xface(i) == centnode
                cons(i,:) = [c_xface(i) 0 0];
            else
                cons(i,:) = [c_xface(i) 0 1];
            end
        end
        for j = 1:numel(n_yface)
            cons(i+j,:) = [n_yface(j) 1 0];
        end
        for k = 1:numel(p_yface)
            cons(i+j+k,:) = [p_yface(k) 1 0];
        end
    end
    
end

%% 3D elements
if type == 2
    
    xmin = min(coords(:,2));
    ymin = min(coords(:,3));
    zmin = min(coords(:,4));
    
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
    
    % 3DOF per node, Ciarletta Shear
     if CCASE == 3
        for i = 1:numel(face_yp)    % top nodes
            cons(i,:)=[face_yp(i) 1 0 0];
        end
        j = size(cons,1);
        for i = 1:numel(face_yn)    % bottom nodes
            cons(j+i,:)=[face_yn(i) 0 0 0];
        end
    end
    
    
    % 3DOF per node, uniaxial tension
    if CCASE == 5
        cons = zeros(numel(face_xc)+numel(face_yn)+numel(face_zn),4);
        for i = 1:numel(face_xc)   % centre nodes
            cons(i,:) = [face_xc(i) 0 1 1];
        end
        for j = 1:numel(face_yn)
            cons(i+j,:) = [face_yn(j) 1 0 1];
        end
        for k = 1:numel(face_zn)
            cons(i+j+k,:) = [face_zn(k) 1 1 0];
        end
    end
    
    if CCASE == 6
        cons = [1 0 0 0;
                2 1 0 0;
                3 0 1 0;
                4 1 1 0;
                5 0 0 1;
                6 1 0 1;
                7 0 1 1];
    end
    
    % Biaxial Tension
    if CCASE == 7
        
        %% Restict planes:
        % set up array full of zeroes
        cons = zeros(numel(face_xc)+numel(face_yc)+numel(face_zc),4);

        % restrict x displacement on central x-face
        for i = 1:numel(face_xc)   % centre nodes
            cons(i,:) = [face_xc(i) 0 1 1];
        end
        % restrict y displacement on central y-face
        for j = 1:numel(face_yc)
            cons(i+j,:) = [face_yc(j) 1 0 1];
        end
        % restrict z displacement on central z-face
        for k = 1:numel(face_zc)
            cons(i+j+k,:) = [face_zc(k) 1 1 0];
        end

        %% symmetric problem:
        % set up array full of zeroes
        cons = zeros(numel(face_xn)+numel(face_yn)+numel(face_zc),4);

        % restrict x displacement on central x-face
        for i = 1:numel(face_xn)   % centre nodes
            cons(i,:) = [face_xn(i) 0 1 1];
        end
        % restrict y displacement on central y-face
        for j = 1:numel(face_yn)
            cons(i+j,:) = [face_yn(j) 1 0 1];
        end
        % restrict z displacement on central z-face
        for k = 1:numel(face_zc)
            cons(i+j+k,:) = [face_zc(k) 1 1 0];
        end

        %% Restrict central node:
        % cons = [centnode 0 0 0];
    end

    if CCASE == 91
       % Constraints for FS mode
       for i = 1:numel(face_xp)    % front face shear in y
            cons(i,:)=[face_xp(i) 0 1 0];
        end
        j = size(cons,1);
        for i = 1:numel(face_xn)    % back face
            cons(j+i,:)=[face_xn(i) 0 0 0];
        end        
    end
        
    if CCASE == 92
        % Constraints for FN mode
       for i = 1:numel(face_xp)    % xp face shear in z
            cons(i,:)=[face_xp(i) 0 0 1];
        end
        j = size(cons,1);
        for i = 1:numel(face_xn)    % xn face
            cons(j+i,:)=[face_xn(i) 0 0 0];
        end        
    end
        
    if CCASE == 93
        % Constraints for SF mode
       for i = 1:numel(face_yp)    % yp face shear in x
            cons(i,:)=[face_yp(i) 1 0 0];
        end
        j = size(cons,1);
        for i = 1:numel(face_yn)    % yn face
            cons(j+i,:)=[face_yn(i) 0 0 0];
        end        
    end
        
    if CCASE == 94
        % Constraints for SN mode
       for i = 1:numel(face_yp)    % yp face shear in z
            cons(i,:)=[face_yp(i) 0 0 1];
        end
        j = size(cons,1);
        for i = 1:numel(face_yn)    % yn face
            cons(j+i,:)=[face_yn(i) 0 0 0];
        end        
    end
        
    if CCASE == 95
        % Constraints for NF mode
       for i = 1:numel(face_zp)    % zp face shear in x
            cons(i,:)=[face_zp(i) 1 0 0];
        end
        j = size(cons,1);
        for i = 1:numel(face_zn)    % zn face
            cons(j+i,:)=[face_zn(i) 0 0 0];
        end        
    end
        
    if CCASE == 96
        % Constraints for NS mode
       for i = 1:numel(face_zp)    % zp face shear in y
            cons(i,:)=[face_zp(i) 0 1 0];
        end
        j = size(cons,1);
        for i = 1:numel(face_zn)    % zn face
            cons(j+i,:)=[face_zn(i) 0 0 0];
        end        
    end  
end

if type == 4

    constraint = 2;

    if constraint == 1
    
        cons = zeros(size(coords));
        x_zero = find(abs(coords(:,2))<=tol);
        y_zero = find(abs(coords(:,3))<=tol);
        z_zero = find(abs(coords(:,4))<=tol);
        % Find all required nodes
        % for i = 1:numel(x_zero)
        %     cons(i,:) = [x_zero(i) 0 1 1];
        % end
        % 
        % for j = 1:numel(y_zero)
        %     cons(j,:) = [y_zero(j) 1 0 1];
        % end
        % 
        % for k = 1:numel(z_zero)
        %     cons(k,:) = [z_zero(k) 1 1 0];
        % end
        % % Define loop to set contrains on each relevant 'face'
        % 
        % 
        % %ignore below, just using to get it running...
        % %cons = [1 0 0 0; 2 0 0 0; 3 0 0 0; 4 0 0 0];
        for i = 1:numel(x_zero)
            cons(x_zero(i),:) = [x_zero(i) 0 1 1];
        end
        
        for j = 1:numel(y_zero)
            if ismember(y_zero(j),x_zero) == 1
                cons(y_zero(j),:) = [y_zero(j) 0 0 1];
            else
                cons(y_zero(j),:) = [y_zero(j) 1 0 1];
            end
        end
        
        for k = 1:numel(z_zero)
            if ismember(z_zero(k),x_zero) == 1
                cons(z_zero(k),:) = [z_zero(k) 0 1 0];
            else
                cons(z_zero(k),:) = [z_zero(k) 1 1 0];
            end
        end

    elseif constraint == 2
        cons = zeros(size(coords));

        nt = mesh_div; %4
        nl = thickness_div + 1; %4
        nm = height_div + 1; %4
        
        %number of nodes in single layer
        nodes_layer = nt * nl; %16
        for i = 1:nt
            i = i + nt*(nl-1);
            for j =1:nm
                cons(i + (j-1)*nodes_layer,:) = [(i + (j-1)*nodes_layer) 0 0 0];
            end
        end
    end
    cons(cons(:,1)==0,:) = [];
        % % find the outside nodes
        % num_theta = 4;
        % num_lambda = 2;
        % 
        % nodes = [1 2 3 4 10 11 12 13];
        % 
        % 
        % for i = 1:numel(nodes)
        %     cons(i,:) = [nodes(i) 0 0 0];
        % end
end