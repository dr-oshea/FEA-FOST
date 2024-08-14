function [DeltaQ] = BFGS_Update(invK0,L,iter,DelQ,R)

% Adjust iter to match the indexing
i = iter;

% Total degrees of freedom
dofT = size(invK0,1);
% L=1:dofT;

% Identity Matrix
dof = size(L,2);
I = eye(dof,dof);

% Calculate b1
R1 = R(:,i);
[v,w] = GetVecs(2,DelQ,R,L);
b1 = (I + v * w.');
for j = 3:i
    [v,w] = GetVecs(j,DelQ,R,L);
    b1 = b1 * (I + v * w.');
end
b1=b1*R1(L);

% Calculate b2
b2 = invK0*b1;

% Calculate new guess
[v,w] = GetVecs(i-0,DelQ,R,L);
result = (I + w * v.');
for j=1:i-2
    [v,w] = GetVecs(i-j,DelQ,R,L);
    result = result * (I + w * v.');
end

DeltaQ = result*b2;

end
    