function CC = Get_CC(F, Lg, Mg, model)

% This script returns the PK2 stress tensor for a given model for use in
% the nonlinear finite element program developed by O'Shea (2018)

% Author:   Daniel O'Shea
% Created:  23 March 2023

% INPUTS:
% F = current deformation gradient
% L = current Lambda tensor
% M = current Mu tensor
% model = marker used to distinguish between various models

% OUTPUTS:
% S = PK2 stress tensor

%% ---------------------------------------------------------------------------
global nset cset
num_param = numel(nset);

C = F.'*F;
I = eye(size(C,1));
IoI = circledot22(I,I);

[V,E] = eig(C);
normE = norm(diag(E));
tol = 10^-3*normE;

% Difference between eigenvalues
d12 = abs(E(1,1)-E(2,2));
d13 = abs(E(1,1)-E(3,3));
d23 = abs(E(2,2)-E(3,3));

% Offset one eigenvalue if two are equal (make them distinct)
if d12<tol
    if E(1,1) >= E(2,2);  E(2,2) = E(2,2) - tol; end
    if E(2,2) >= E(1,1);  E(1,1) = E(1,1) - tol; end
end
if d13<tol
    if E(1,1) >= E(3,3);  E(3,3) = E(3,3) - tol; end
    if E(3,3) >= E(1,1);  E(1,1) = E(1,1) - tol; end
end
if d23<tol
    if E(3,3) >= E(2,2);  E(2,2) = E(2,2) - tol; end
    if E(2,2) >= E(3,3);  E(3,3) = E(3,3) - tol; end
end
Ct = V*E/V;
RR = circledot22(Ct,I)-circledot22(I,Ct)+circlecross22(Ct,I)-circlecross22(I,Ct)+circlecross22(Ct*Ct,Ct*Ct);


% Lambda Term:
invC = inv(C);
lnC = logm(C);
E0 = (1/2)*lnC;
invC_C = ITF(invC,-invC*invC,C);
E0_C = ITF(E0,(1/2)*invC,C);
lnC_C = ITF(log(C),invC,C);
G = DblCon42(Lg,E0);
G_C = DblCon44(Lg,E0_C);
BigTerm = SinCon42(lnC_C,G)-SinCon24(G,lnC_C) + SinCon24(lnC,G_C) - SinCon42(G_C,lnC) - G_C + SinCon42(invC_C,G*C) + SinCon24(invC,SinCon42(G_C,C)) + SinCon24(invC,SinCon24(G,IoI));
CC_lam = (SinCon42(invC_C,G)+SinCon24(invC,G_C)) + DblCon44(Inv4(Transpose4(RR,'dot'),'dot'),BigTerm);
Sbar_L = DblCon24(G,2*E0_C-circledot22(invC,I));
SL = 2*DblCon24(G,E0_C);

%CC_lam = DblCon44(circledot22(Power2(C,-1),I),Lg,'dot');


CC_mu = zeros(3,3,3,3); Sbar_M = zeros(3,3); SM = Sbar_M;
for z = 1:num_param
    ni = nset(z);
    ci = cset(z);
    En = (1/ni)*(Power2(C,ni/2)-I);
    Cn2 = Power2(C,ni/2);
    Cn21 = Power2(C,ni/2-1);
    Cn2d = (ni/2)*Cn21;
    Cn21d = (ni/2-1)*Power2(C,ni/2-2);
    G = DblCon42(ci*Mg,En);

    En_C = ITF(En,(1/2)*Cn21,C);
    Cn2_C = ITF(Cn2,Cn2d,C);
    Cn21_C = ITF(Cn21,Cn21d,C);
    G_C = DblCon44(ci*Mg,En_C);

    BigTerm = (1-ni/2)*(SinCon42(Cn2_C,G)+SinCon24(Cn2,G_C))- SinCon24(G,Cn2_C) - SinCon42(G_C,Cn2) + (ni/2)*(SinCon42(Cn21_C,G*C) + SinCon24(Cn21,SinCon42(G_C,C)) + SinCon24(Cn21,SinCon24(G,IoI)));
    CC_mu = CC_mu + (SinCon42(Cn21_C,G)+SinCon24(Cn21,G_C)) + (2/ni)*DblCon44(Inv4(Transpose4(RR,'dot'),'dot'),BigTerm);

    %CC_mu = CC_mu + DblCon44(circledot22(Power2(C,ni/2-1),I),ci*Mg,'dot');

    Sbar_M = Sbar_M + DblCon24(G,2*En_C-circledot22(Cn21,I));
    SM = SM + 2*DblCon24(G,En_C);
end
Sbar = Sbar_L + Sbar_M;
Bonus = circlestar22(I,Sbar.')-circlestar22(Sbar,I)-trace(Sbar)*IoI+trace(C*C*Sbar)*(circledot22(C,I)+circledot22(I,C))+circlecross22(C*C,C*Sbar+Sbar*C);
Bonus_ = -DblCon44(Inv4(Transpose4(RR,'dot'),'dot'),Bonus);

CC = 2*(CC_lam + CC_mu + Bonus_);     % factor of 2 because CC = 2*(dS/dC);

S = SL + SM;

end