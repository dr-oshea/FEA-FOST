function Q = SolveQ(K,P,L,input)

% This script solves the displacement vector by inverting the stiffness
% matrix and multiplying by the load vector

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

nodes = size(input.ND,1);
dofN = size(input.ND,2)-1;

dofT = nodes*dofN;

%%

Q = zeros(1,dofT);

% Calculate displacements of unconstrained degrees of freedoms
Q_ = K(L,L)\P(L)';

% Form complete displacement vector 
Q(L) = Q_;

end

