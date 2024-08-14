function [dofN, ndsE, strC] = input_info(type)

% This function returns element information depending on the element
% selected for use

% Author:   Daniel O'Shea
% Created:  12 April 2018

% INPUTS:
% type = element type (1 = 2D quad, 2 = 3D hex)

% OUTPUTS:
% dofN = degrees of freedom per node
% ndsE = nodes per element
% strC = stress/strain components

%% ------------------------------------------------------------------------

% 2D Quadrilateral Element
if type == 1 || type == 3 
    
    dofN = 2;
    ndsE = 4;
    strC = 4;
    
% 3D Hexahedral Element
elseif type == 2 || type == 4
    
    dofN = 3;
    ndsE = 8;
    strC = 9;

% IFT 2D Brick Element
elseif type == 3
    
    % new element
    
    
    
end




end