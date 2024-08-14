function [S,CC] = Get_S_CC(F, modelData, fibres)

constitutiveModel = modelData(1);   % eventually carry this from problem definition...

% Preliminary kinematic tensors
dimensions = size(F,1);
I = eye(dimensions);

J = det(F);         % Jacobian (must be greater than or equal to 1)
C = F.'*F;      % right Cauchy Green tensor
E2 = (1/2)*(C-I);   % Green-Lagrange strain tensor

% if J < 1+10^(-5) && J > 1-10^(-5)
%     lnJ = 0;
% else
%     lnJ = log(J);
% end
lnJ = log(J);
invC = inv(C);

switch constitutiveModel
    
    case 1    % Myocardium model (published 2022)
        
        % retrive testing parameters
        L00 = modelData(2);
        LFF = modelData(3); 
        LSS = modelData(4); 
        LNN = modelData(5);
        M00 = modelData(6);
        MFS = modelData(7);
        kappa1 = modelData(8); 
        kappa2 = modelData(9);
        k = modelData(10);
        if dimensions == 2
            LNN = 0;    % remove normal fibre for 2D problem.
        end
        
        % dont alter anything below here.
        
        % structural tensors
        if isempty(fibres)
            M1 = zeros(dimensions,dimensions); M2 = M1; M3 = M1;
            M1(1,1) = 1;
            M2(2,2) = 1;
            if dimensions == 3
                M3(3,3) = 1;
            end
        else
            mf = fibres.f';
            ms = fibres.s';
            mn = fibres.n';
            M1 = circlecross11(mf,mf);
            M2 = circlecross11(ms,ms);
            M3 = circlecross11(mn,mn);
        end

        % Second order structural tensors
        M1 = kappa1 * I + (1 - 3*kappa1) * M1;
        M2 = kappa2 * I + (1 - 3*kappa2) * M2;
        
        % Lambda Tensor
        LT = LFF*circlecross22(M1,M1)+LSS*circlecross22(M2,M2)+LNN*circlecross22(M3,M3);
        
        % PK2 Stress Tensor
        S = L00*lnJ*invC + M00*E2 + (k-1)*trace(E2*E2)^(k-2)*(LFF*trace(M1*E2)^2 + LSS*trace(M2*E2)^2 + LNN*trace(M3*E2)^2)*E2 + (trace(E2*E2)^(k-1))*(LFF*trace(M1*E2)*M1 + LSS*trace(M2*E2)*M2 + LNN*trace(M3*E2)*M3);
        % S = M00*E2 + (k-1)*trace(E2*E2)^(k-2)*(LFF*trace(M1*E2)^2 + LSS*trace(M2*E2)^2 + LNN*trace(M3*E2)^2)*E2 + (trace(E2*E2)^(k-1))*(LFF*trace(M1*E2)*M1 + LSS*trace(M2*E2)*M2 + LNN*trace(M3*E2)*M3);
        

        % Tangent Elasticity Tensor
        if F == I
            CC = L00*circlecross22(I,I) + M00*circledot22(I,I);
        else
            CC = L00*(circlecross22(invC,invC)-2*lnJ*circledot22(invC,invC)) + M00*circledot22(I,I) ...
            + trace(E2*E2)^(k-1)*QuadCon44(LT,circlecross22(E2,E2))*(k-1)*trace(E2*E2)^(k-3)*((k-2)*circlecross22(E2,E2) + (1/2)*trace(E2*E2)*circledot22(I,I)) ...
            + (k-1)*trace(E2*E2)^(k-2)*(circlecross22(E2,DblCon42(LT,E2))+circlecross22(DblCon42(LT,E2),E2)) ...
            + trace(E2*E2)^(k-1)*circledot22(I,I);
        end
        
    case 2      % Gultekin model

        % retrive testing parameters
        a = modelData(2);
        b = modelData(3); 
        af = modelData(4); 
        bf = modelData(5);
        as = modelData(6);
        bs = modelData(7);
        afs = modelData(8);
        bfs = modelData(9);
        kappaF = modelData(10); 
        kappaS = modelData(11);
        
        % structural tensors (hard-coded)        
        if isempty(fibres)
            M1 = zeros(dimensions,dimensions); M2 = M1; M3 = M1;
            M1(1,1) = 1;
            M2(2,2) = 1;
            if dimensions == 3
                M3(3,3) = 1;
            end
        else
            mf = fibres.f';
            ms = fibres.s';
            mn = fibres.n';
            M1 = circlecross11(mf,mf);
            M2 = circlecross11(ms,ms);
            M3 = circlecross11(mn,mn);
        end
        
        % Second order structural tensors
        H1 = kappaF * I + (1 - 3*kappaF) * M1;
        H2 = kappaS * I + (1 - 3*kappaS) * M2;
        HFS = (1/2)*circlecross11([1;0;0],[0;1;0]) + circlecross11([1;0;0],[0;1;0]);
        
        % Fourth-order structural tensors
        HHF = circlecross22(H1, H1);
        HHS = circlecross22(H2, H2);
        HHFS = circlecross22(HFS, HFS);
        
        % PK2 Stress Tensor
        I1 = trace(C*I);
        I4f = trace(C*H1);
        I4s = trace(C*H2);
        I8fs = trace(C*HFS);
        S = 2*(a/2*exp(b*(I1 - 3))*I + af*exp(bf*(I4f - 1)^2)*(I4f - 1)*H1 + af*exp(bf*(I4s - 1)^2)*(I4s - 1)*H2 + afs*exp(bfs*(I8fs^2))*I8fs*HFS);
            
        % Tangent Elasticity Tensor
        if F == I
            psi1 = a*b/2;
            psi4 = af;
            psi6 = as;
            psi8 = afs;
        else
            psi1 = a*b/2*exp(b*(I1 - 3));
            psi4 = af*(1 + 2*bf*(I4f - 1) ^2)*exp(bf*(I4f - 1)^2);
            psi6 = as*(1 + 2*bs*(I4s - 1) ^2)*exp(bs*(I4s - 1)^2);
            psi8 = afs*(1 + 2*bfs*I8fs^2)*exp(bfs*I8fs^2);
        end
        CC = 4*(psi1*circlecross22(I,I) + psi4*HHF + psi6*HHS + psi8*HHFS);
    

    case 3      % Isotropic SVK model

        % retrive testing parameters
        E = modelData(2);
        v = modelData(3); 
        lambda = E*v/((1+v)*(1-2*v));
        mu = E/(1+v);

        S = lambda * trace(E2) * I + mu * E2;

        CC = lambda * circlecross22(I,I) + mu * circledot22(I,I);


    otherwise
        
        fprintf('Unable to calculate S or CC for Tangent Stiffness. Model is not defined');
        
end


end