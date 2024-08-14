function PN = Final_Load_Vector(input)

% This script converts the array of prescribed nodal loads in to the
% "final" load vector for the nonlinear elasticity problem

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% input = structure containing all input data

% OUTPUTS:
% PN = final load vector for problem

%% -------------------------------------------------------------------------

nodes = size(input.ND,1);
dofPerNode = size(input.ND,2)-1;

totalDOFS = nodes*dofPerNode;


PN = spalloc(1,totalDOFS,1);

for i=1:size(input.LOAD,1)
    
    t4 = input.LOAD(i,1)*dofPerNode-(dofPerNode-1);
    
    for k=1:size(input.LOAD,2)-1
    
        PN(t4+k-1) = input.LOAD(i,k+1);
    
    end
    
end
