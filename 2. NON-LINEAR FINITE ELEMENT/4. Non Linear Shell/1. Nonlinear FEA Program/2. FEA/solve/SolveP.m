function P = SolveP(Q,input,ndsPerElem)

% This script solves the load vector by means of a integration a nonlinear
% function based on the orthotropic SEF proposed by O'Shea et. al. (2018)

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% K = stiffness matrix
% P = load vector
% L = array containing unconstrained degrees of freedom
% input = contains all input information for the problem

% OUTPUT:
% Q = displacement vector

%% ------------------------------------------------------------------------

dofN = size(input.ND,2)-1;
type = input.EL(1,2);

%%

P = Assemble_Load_Vector(Q,input,dofN,ndsPerElem,type);

end

