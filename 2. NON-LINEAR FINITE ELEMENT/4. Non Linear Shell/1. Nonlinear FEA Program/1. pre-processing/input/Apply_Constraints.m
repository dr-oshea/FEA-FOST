function L = Apply_Constraints(input)

% This script constrains the global system

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% input = contains all input information for problem

% OUTPUTS:
% L = vector containing


%% -------------------------------------------------------------------------

nodes = size(input.ND,1);
dofN = size(input.ND,2)-1;

dofT = nodes*dofN;

% Initially, L is list of all DOFs
L0 = 1:dofT;

ind = 1;

% For every constrained node
for i=1:size(input.CON,1)
    
    % For each degree of freedom
    for k=1:dofN
        
        % If the DOF is constrained
        if input.CON(i,1+k)==0
            
            % Store DOF in terms of global DOF list
            L(ind) = input.CON(i,1)*dofN-(dofN-k);
            ind = ind+1;
            
        end
    end
end

% Removes DOFs that are constrained from list of global DOFs
L = setdiff(L0,L);  