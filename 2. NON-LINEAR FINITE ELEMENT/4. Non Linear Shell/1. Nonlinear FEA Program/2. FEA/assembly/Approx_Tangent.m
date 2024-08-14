function KT = Approx_Tangent(PN,PO,QN,QO,L,dofT)

% This function uses the finite difference method to approximate the
% tangent stiffness based on the last two known displacement solutions

% Author:   Daniel O'Shea
% Created:  11 April 2018

% INPUTS:
% KN = secant stiffness at current point
% KO = secant stiffness at previous point
% dofT = total degrees of freedom

% OUTPUTS:
% KT = two-point approximation of tangent stiffness

%% ------------------------------------------------------------------------

num = numel(L);

KT = zeros(dofT,dofT);

for ii = 1:dofT
    i = ii;
    
    for jj = 1:num
        j = L(jj);
        
        KT(i,j) = KT(i,j) + (PN(i) - PO(i)) / (QN(j) - QO(j));
        
    end
    
end

KT = sparse(KT);