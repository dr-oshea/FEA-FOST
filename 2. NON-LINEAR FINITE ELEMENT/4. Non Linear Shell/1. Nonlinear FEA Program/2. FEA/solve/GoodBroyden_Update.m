function [DeltaQ] = GoodBroyden_Update(invK,R1,Q1,Pn)

% Give a dummy updated estimate using current secant stiffness
DeltaQ2_ = invK*R1(L)';
DeltaQ2 = zeros(1,dofT);
DeltaQ2(L) = DeltaQ2_;
Q2 = Q1 + DeltaQ2';
P2 = SolveP(Q2,input,ndsE);
R2 = full(Pn)-P2';

% Vectors required
yy = (R2(L)-R1(L))';
ss = (Q2(L)-Q1(L));

% Update inverse secant stiffness
invK = invK + (ss-invK*yy)*(ss.'*invK)/(ss.'*invK*yy);
invK = -invK;

% Improved estimate and residual
DeltaQ = invK*R1';
        
end
    