function coords = input_nodal(mesh_div,aspect,xdim,ydim,zdim,type,thickness_div,height_div)

% This function generates an array of nodal coordinate information for a
% specified two-dimensional geometry. Currently limited to rectangular
% meshing

% Author:   Daniel O'Shea
% Created:  16 March 2018

% INPUTS:
% mesh_div = nodal divisions along x-dimension
% aspect = aspect ratio determining divisions along y-dimension
% xdim = total magnitude of geometry along x-dimension
% ydim = total magnitude of geometry along y-dimension
% zdim = total magnitude of geometry along z-dimension
% type = element type

% OUTPUTS:
% coords = array storing nodal information for discretised mesh

% UPDATE LOG:
% 12/04/2018 - included support for a 3D element



%% ---------------------------------------------------------------------------
% Determine element size


% For a CSTQ element
if type == 1  
    spacing = xdim/(mesh_div);         
    xspace = spacing;
    yspace = aspect*spacing;
    
    % Number of nodes along each dimension
    nx = xdim/xspace+1;
    ny = ydim/yspace+1;

    % Total nodes in discretisation
    totalnodes = nx*ny;

    % Determine coordinates of each node
    coords=zeros(totalnodes,3);

    for i=1:totalnodes
            x=rem(i+nx-1,nx)*xspace;
            y=rem(ceil(i/nx)-1,ny)*yspace;
            coords(i,:)=[i x y];
    end

end

% For a CSTQ element in polar coordinates
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
    lambdaMin = zdim;
    lambdaMax = zdim + ydim;
    
    % define arrays of points in each polar coordinate
    % divTheta = mesh_div + 1;
    % divRadius = mesh_div * aspect;
    divTheta = mesh_div + 1;
    divLambda = mesh_div * aspect;

    theta = linspace(thetaMin, thetaMax, divTheta);
    lambda = linspace(lambdaMin, lambdaMax, divLambda);

    % Number of nodes along each dimension
    nt = numel(theta);
    nl = numel(lambda);

    % Total nodes in discretisation
    totalnodes = (nt-1)*nl;

    % Initialise array for coordinates of each node
    coords=zeros(totalnodes,3);

    % Loop throough each polar coordinate and calculate Cartesian coordinate
    a = 1;
    for i=1:nl 
        for j=1:(nt-1)  % leave out last theta as is duplicate node
            t = theta(j);
            l = lambda(i);
            x = l*cos(t);
            y = l*sin(t);
            coords(a,:) = [a, x, y];
            a = a + 1;
        end
    end

end

% For a BRICK element
if type == 2
    spacing = zdim/(mesh_div);         
    xspace = aspect*spacing;
    yspace = aspect*spacing;
    zspace = spacing;
    
    % Number of nodes in each direction
    nx = xdim/xspace+1;
    ny = ydim/yspace+1;
    nz = zdim/zspace+1;
    
    % Total nodes in discretisation
    totalnodes = nx*ny*nz;
    
    % Determine coordinates of each node
    coords=zeros(totalnodes,4);
    for i=1:totalnodes
            x=rem(i+nx-1,nx)*xspace;
            y=rem(ceil(i/nx)-1,ny)*yspace;
            z=(ceil(i/(nx*ny)-1))*zspace;
            coords(i,:)=[i x y z];
    end
    
end

% For a CSTQ element in polar coordinates
if type == 4

    thetaMin = 0;
    thetaMax = xdim;
    lambdaMin = zdim;
    lambdaMax = zdim + ydim;
    muMin = pi/3;  %0.01 % pi/3
    muMax = pi/2;
    d = 60;
    
    % define arrays of points in each polar coordinate
    % divTheta = mesh_div + 1;
    % divRadius = mesh_div * aspect;
    divTheta = mesh_div + 1;
    divLambda = thickness_div + 1;
    divMu = height_div + 1;

    theta = linspace(thetaMin, thetaMax, divTheta);
    lambda = linspace(lambdaMin, lambdaMax, divLambda);
    mu = linspace(muMax, muMin, divMu); %intentionally put backwards to generate nodes starting from open end
    % mu = exp(linspace(log(muMax),log(muMin),divMu)); %logarithmically spaced mu values (more towards base of ventricle)
    linspace(muMax, muMin, divMu); %intentionally put backwards to generate nodes starting from open end

 
    % Number of nodes along each dimension
    nt = numel(theta);
    nl = numel(lambda);
    nm = numel(mu);

    % Total nodes in discretisation
    totalnodes = (nt-1)*nl*nm;

    % Initialise array for coordinates of each node
    coords=zeros(totalnodes,4);

    % Loop throough each polar coordinate and calculate Cartesian coordinate
    a = 1;
    for i=1:nm 
        for j=1:nl  % leave out last theta as is duplicate node
            for k=1:(nt-1)
                m = mu(i);
                l = lambda(j);
                t = theta(k);  %right now its labelling nodes down the heart, and then around theta and then going outwards. need to make it go around theta, then outwards and then updwards

                x = d*cosh(l)*cos(m);
                y = d*sinh(l)*sin(m)*cos(t);
                z = d*sinh(l)*sin(m)*sin(t);
                coords(a,:) = [a, x, y, z];
                a = a + 1;
            end
        end
    end

end
    
end