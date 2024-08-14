function B = B_Matrix_2D(dN)

% This script returns the B1 matrix for a 2D quad element. this defines the
% relationship between displacements and displacement gradients

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% dN = array storing shape function derivatives in global coords

% OUTPUTS:
% B1 = B1 matrix relating displacement gradients to displacements

%% ------------------------------------------------------------------------

% CSTQ Element
ndsPerElem = 4;
dofN = 2;
strC = 4;

for n = 1:ndsPerElem
    
    row0 = 1;
    rowf = strC;
    col0 = dofN*(n-1)+1;
    colf = dofN*(n-1)+dofN;
    
    B(row0:rowf,col0:colf)= [    dN(1,n)        0    ; 
                                    0        dN(2,n) ;
                                 dN(2,n)        0    ;
                                     0        dN(1,n)];

end