function B = B_Matrix_3D(dN)

% This script returns the B1 matrix for a 3D brick element. This defines 
% the relationship between displacements and displacement gradients

% Author:   Daniel O'Shea
% Created:  12 April 2018

% INPUTS:
% dN = array storing shape function derivatives in global coords

% OUTPUTS:
% B1 = B1 matrix relating displacement gradients to displacements

%% ------------------------------------------------------------------------

% Brick Element
ndsPerElem = 8;
dofN = 3;
strC = 9;
    
for n = 1:ndsPerElem
    
    row0 = 1;
    rowf = strC;
    col0 = dofN*(n-1)+1;
    colf = dofN*(n-1)+dofN;
    
    B(row0:rowf,col0:colf)= [    dN(1,n)        0           0   ; 
                                    0        dN(2,n)        0   ;
                                    0           0        dN(3,n);
                                    0        dN(3,n)        0   ;
                                    0           0        dN(1,n);
                                 dN(2,n)        0           0   ;
                                    0           0        dN(2,n);
                                 dN(3,n)        0           0   ;
                                    0        dN(1,n)        0   ];

end