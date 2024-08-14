%% Equations for a 2D plate element
dim = 2;

syms u1 u2

% We know the values of u1 and u2 by interpolating with values at all 4
% nodes of the element

syms u1_1 u1_2 u2_1 u2_2

% We know the values of derivatives of displacement components by using the
% derivative of shape functions and interpolating

%% CAUCHY-GREEN DEFORMATION TENSOR
delta = [1 0; 0 1];

gradu = [u1_1 u1_2; u2_1 u2_2];

C = sym(zeros(dim,dim));
for i = 1:dim
    for j = 1:dim
        A = gradu(i,j) + gradu(j,i) + delta(i,j);
        B = 0;
        for k = 1:dim
            B = B + gradu(k,i)*gradu(k,j);
        end
        C(i,j) = A + B;
    end
end

%% EIGENDECOMPOSITION

[V, D] = eig(C);