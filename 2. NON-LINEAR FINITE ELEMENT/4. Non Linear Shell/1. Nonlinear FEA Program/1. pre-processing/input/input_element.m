function [element, fibres] = input_element(coords,mesh_div,aspect,layers,xdim,ydim,zdim,theta1,theta2,type,thickness_div,height_div,use_dispersion)

%%  input_element

% This script generates an array containing element connectivity
% information using paramters defined in FEA_nonlin_???.m input files. The 
% elements are restricted to square (2D)

% Author:   Daniel O'Shea
% Created:  16 March 2018

% INPUTS:
% coords = nodal coordinate data
% mesh_div = element divisions along x-dimension of element (2D models)
% layers = number of layers for a composite plate
% xdim = dimension of element in the x-direction [mm]
% ydim = dimension of element in the y-direction [mm]
% theta1 = defined fibre orientation #1
% theta2 = defined fibre orientation #2

% OUTPUTS:
% element = contains all appropriate data related to all local elements

% UPDATE LOG:
% 12/04/2018 - included support for a 3D element
% 13/10/2023 - included method to define FSN orientation vectors 



%% ---------------------------------------------------------------------------
tol = 10^-6;
fibres = [];
fvectors = []; svectors = []; nvectors = [];

if type == 1
    
    % Determine element size
    spacing = xdim/(mesh_div);         
    xspace = spacing;   % alter aspect ratio of elements
    yspace = spacing;

    % Number of nodes in each direction
    x = xdim/xspace+1;
    y = ydim/yspace+1;
    
    % Total number of elements
    totalelements = (x-1)*(y-1);
    xele = x-1;
    yele = y-1;
    
    xmat=zeros(1,xele); ymat=zeros(1,yele);
    for i=1:xele
        xmat(:,i)=(i-1)*xspace; end
    for i=1:yele
        ymat(:,i)=(i-1)*yspace; end

    element = zeros(totalelements,7);
    elementmatrix = zeros(totalelements,3);
    for j=1:totalelements           
        xval = xmat(:,1+rem(j+xele-1,xele));
        yval = ymat(:,1+rem(ceil(j/xele)-1,yele));
        elementmatrix(j,:) = [j xval yval];  
        node1 = find(coords(:,2)==elementmatrix(j,2)&coords(:,3)==elementmatrix(j,3));
        node2 = node1+1;
        node3 = node1+x+1;
        node4 = node3-1;

        if j<=(1/layers)*totalelements
            theta = theta1; 
        elseif j<=(2/layers)*totalelements
            theta = theta2;
        end

        element(j,:)=[j type node1 node2 node3 node4 theta];

    end
end

if type == 3

    % The program defines 'xdim', 'ydim', 'zdim', in the RunFEA file to
    % define geometry, so easiest to continue to give them meaning here but 
    % related to the polar coordinates.

    % Let xdim = maximum theta values = 2*pi (complete circle)
    % Let ydim = thickness of annulus = outer radius - inner radius
    % Let zdim = inner radius
    
    % Limits of polar coordinates
    thetaMin = 0;
    thetaMax = xdim;
    radiusMin = zdim;
    radiusMax = zdim + ydim;
    
    % define arrays of points in each polar coordinate
    theta = linspace(thetaMin, thetaMax, mesh_div + 1);
    radius = linspace(radiusMin, radiusMax, mesh_div);

    % Number of nodes along each dimension
    nt = numel(theta);
    nl = numel(radius);     
    
    % Total number of elements
    numElementsT = nt-1;
    numElementsL = nl-1;
    totalElements = numElementsT * numElementsL;
    
    for a = 1:totalElements
        node1 = a;  % based on node numbering convention, 'node 1' number should always be same as element number
        if mod(node1,(nt-1)- 0) == 0  % account for final element each cycle of theta
            node2 = node1 - ((nt-1)-1) + 0;
        else
            node2 = node1 + 1;      
        end
        node3 = node2 + (nt-1) - 0;     % shifted 'outwards' one cycle of theta  %changed from "(nt-1) - 1" to "(nt-1) - 0"
        node4 = node1 + (nt-1) - 0;     % shifted 'outwards' one cycle of theta
        element(a,:) = [a, type, node1, node2, node3, node4, theta1];
    end

end


%% 3D element
if type == 2
    
    % Determine element size
    spacing = zdim/(mesh_div);         
    xspace = aspect*spacing;   % alter aspect ratio of elements (aspect: 1 = cubic)
    yspace = aspect*spacing;
    zspace = spacing;

    % Number of nodes in each direction
    x = xdim/xspace+1;
    y = ydim/yspace+1;
    z = zdim/zspace+1;
    
    % Total number of elements in mesh
    totalelements = (x-1)*(y-1)*(z-1);
    xele = x-1;
    yele = y-1;
    zele = z-1;

    xmat=zeros(1,xele); ymat=zeros(1,yele); zmat=zeros(1,zele);
    for i=1:xele
        xmat(:,i)=(i-1)*xspace; end
    for i=1:yele
        ymat(:,i)=(i-1)*yspace; end
    for i=1:zele
        zmat(:,i)=(i-1)*zspace; end

    element=zeros(totalelements,11);
    elementmatrix=zeros(totalelements,4);
    for j=1:totalelements
        xval=xmat(:,1+rem(j+xele-1,xele));
        yval=ymat(:,1+rem(ceil(j/xele)-1,yele));
        zval=zmat(:,ceil(j/(xele*yele)));
        elementmatrix(j,:)=[j xval yval zval];  
        
        % Relate local node numbers to global node numbers
        node1 = find(abs(coords(:,2)-elementmatrix(j,2))<tol&abs(coords(:,3)-elementmatrix(j,3))<tol&abs(coords(:,4)-elementmatrix(j,4))<tol);
        node2 = node1+1;
        node3 = node1+x+1;
        node4 = node3-1;
        node5 = node1+x*y;
        node6 = node5+1;
        node7 = node5+x+1;
        node8 = node7-1;

        if j<=(1/layers)*totalelements
            theta = theta1; 
        elseif j<=(2/layers)*totalelements
            theta = theta2;
        end

        element(j,:)=[j type node1 node2 node3 node4 node5 node6 node7 node8 theta];

    end

end

if type == 4       % 3D ellipsoid
    
    % 2D code to be changed to 3d code....


    % The program defines 'xdim', 'ydim', 'zdim', in the RunFEA file to
    % define geometry, so easiest to continue to give them meaning here but 
    % related to the polar coordinates.

    % Let xdim = maximum theta values = 2*pi (complete circle)
    % Let ydim = thickness of annulus = outer radius - inner radius
    % Let zdim = inner radius
    
    % Limits of polar coordinates
    thetaMin = 0;
    thetaMax = xdim;
    lambdaMin = zdim;
    lambdaMax = zdim + ydim;
    muMin = 0.01;
    muMax = pi/2;
    
    divTheta = mesh_div + 1;
    divLambda = thickness_div + 1;
    divMu = height_div + 1;

    % define arrays of points in each polar coordinate
    theta = linspace(thetaMin, thetaMax, divTheta);
    lambda = linspace(lambdaMin, lambdaMax, divLambda);
    mu = linspace(muMin,muMax,divMu);

    % Number of nodes along each dimension
    nt = numel(theta);
    nl = numel(lambda);
    nm = numel(mu);
    
    % Total number of elements
    numElementsT = nt-1; %3
    numElementsL = nl-1; %2
    numElementsM = nm-1; %2
    ElementsPerLayer = numElementsT*numElementsL; %6
    totalElements = numElementsT * numElementsL * numElementsM;
    
    layer = 0;     % tells you which value of 'mu' you are one (i.e one full ring)
    layer_lambda = 0;  % tells you which value of lambda you are on for some given mu
    alpha_nodes = linspace(0,180,nl); %0 to 180
    %alpha_nodes = linspace(90,90,n1);
    alpha_elements = mean([alpha_nodes(1:end-1); alpha_nodes(2:end)]);
    fvectors = zeros(totalElements,3);
    svectors = zeros(totalElements,3);
    nvectors = zeros(totalElements,3);
    for a = 1:totalElements
        if mod(a-1,numElementsT) == 0
            layer_lambda = layer_lambda + 1;
        end
        if mod(a-1,ElementsPerLayer) == 0
            layer = layer + 1;
            layer_lambda = 1;
        end
        node1 = a + (layer - 1)*numElementsT;  % based on node numbering convention, 'node 1' number should always be same as element number
        if mod(node1, (nt-1)) == 0  % account for final element each cycle of theta
            node2 = node1 - ((nt-1)-1);
        else
            node2 = node1 + 1;      
        end
        node3 = node2 + (nt-1);     % shifted 'outwards' one cycle of theta
        node4 = node1 + (nt-1);     % shifted 'outwards' one cycle of theta
        node5 = node1 + (nt - 1)*nl;
        node6 = node2 + (nt - 1)*nl;
        node7 = node3 + (nt - 1)*nl;
        node8 = node4 + (nt - 1)*nl;

        % node5-node8 will be node1-node4 shifted by the number of elements
        % in a single value of mu (I think..) might be adding nr*(nt-1),
        % i.e. the number of elements in one donut

        element(a,:) = [a, type, node1, node2, node3, node4, node5, node6, node7, node8, theta1];
        % just leave 'theta1' as last entry so that the code works (it
        % doesnt do anything currently)
        
        % Change below to have x, y, z components of vectors
        
        
        %Getting sheet direction (I THINK THE WHOLE TIME I HAD SHEET AND
        %SHEET NORMAL IN THE WRONG DIRECTIONS!!!!!!)

        % sheet = cross((coords(node6,2:4) - coords(node5,2:4)),(coords(node1,2:4) - coords(node5,2:4)));
        % sheet_unit = sheet/norm(sheet);
        sheetnorm = cross((coords(node6,2:4) - coords(node5,2:4)),(coords(node1,2:4) - coords(node5,2:4)));
        sheetnorm_unit = sheetnorm/norm(sheetnorm);
        
        if use_dispersion == 1
            alpha = 90;
        else
            alpha = alpha_elements(layer_lambda);
        end

        fibre_vertical = cross(sheetnorm,(coords(node6,2:4) - coords(node5,2:4)));
        fibre_vertical_unit = fibre_vertical/norm(fibre_vertical);
        fibre = cosd(alpha)*fibre_vertical_unit + cross(sheetnorm_unit,fibre_vertical_unit)*sind(alpha); %result of this equation is a unit vector

        % sheetnorm = cross(sheet,fibre);
        % sheetnorm_unit = sheetnorm/norm(sheetnorm);

        sheet = cross(sheetnorm,fibre);
        sheet_unit = sheet/norm(sheet);

        % alpha_values = linspace(...)
        % alpha = alpha_values(layer_theta)

        % Debugging from Daniel
        %fibre = [1, 0, 0];
        %sheet_unit = [0, 1, 0];
        %sheetnorm_unit = [0, 0, 1];


        fvectors(a,:) = fibre;
        svectors(a,:) = sheet_unit; %[0, 1, 0];
        nvectors(a,:) = sheetnorm_unit; %[0, 0, 1];
        
    end

end

if ~isempty(fvectors) && ~isempty(svectors) && ~isempty(nvectors)
    fibres.f = fvectors;
    fibres.s = svectors;
    fibres.n = nvectors;
end