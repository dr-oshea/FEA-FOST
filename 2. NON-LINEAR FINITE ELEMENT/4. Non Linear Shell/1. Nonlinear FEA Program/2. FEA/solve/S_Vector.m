function [S,E0v,Env] = S_Vector(Cvec,L,M)

% This script determines the second Piola-Kirchhoff stress tensor for a
% given displacement field and constitutive relationship. It is returned in
% a vectorised form. the formulations present here are based on those
% presented in the O'Shea, Attard, Kellermann paper (2018)

% Author:   Daniel O'Shea
% Created:  21 March 2018

% INPUTS:
% q = local displacement vector
% C = right Cauchy-Green deformation tensor (vectorised)
% L = lambda material matrix
% M = mu material matrix

% OUTPUTS:
% S = stress vector

% UPDATE LOG:
% 13/04/2018 - Added funcitonality for 3D problems
% 26/04/2018 - Added strain tensors as output

%% ------------------------------------------------------------------------
global nset cset
num_param = numel(nset);

% Dimensions
dim = size(Cvec,1)^0.5;

% Identity Tensor
I = Iden2(dim);

% Right Cauchy-Green deformation tensor
C = V2M(Cvec);

% Logarithmic Strain and Derivative
E0 = (1/2)*Func2(C,'log');
E0v = M2V(E0);
E0_C = (1/2)*ITF(logm(C),inv(C),C);

% Polynomial Strain and Derivative
En = zeros(dim,dim,num_param);
Env = zeros(dim*dim,num_param);
En_C = zeros(dim,dim,dim,dim,num_param);
for i = 1:num_param
    ni = nset(i);
    En(:,:,i) = (1/ni)*(Power2(C,ni/2)-I);
    Env(:,i) = M2V(En(:,:,i));
    En_C(:,:,:,:,i) = (1/ni)*ITF(Power2(C,ni/2),(ni/2)*Power2(C,ni/2-1),C);
end

% Lambda Stresses:
CL = 2 * T2M(E0_C,'dot') * L ;
SL = CL * E0v;

% Mu Stresses
SM = zeros(dim*dim,1);
for i = 1:num_param
    ci = cset(i);
    CM = 2 * T2M(En_C(:,:,:,:,i),'dot') * (ci .* M) ;
    SM = SM + (CM * Env(:,i));
end

% Stress Vector:
S = SL + SM;

end