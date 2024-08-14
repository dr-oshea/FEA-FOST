function Ke = Element_Tangent_Stiffness_2D (q,elCoords, theta, L, M)

% This script generates the element tangent stiffness matrix at a 
% given displacement field for a 2D quad element. The integral is 
% determined using Gaussian integration

% Author:   Daniel O'Shea
% Created:  22 March 2022

% INPUTS:
% q = local vector of known displacements
% elCoords = nodal cordinates for local nodes
% theta = material orientation

% OUTPUTS:
% K0e = element stiffness matrix

%% ---------------------------------------------------------------------------
ndsPerElem = 4;
dofN = 2;
phi = deg2rad(theta);
I = Iden2(2);
IoI = circledot22(I,I);

% Global Element Constitutive Law
[T4,~,~,~,~] = Transform2(phi);
Lg = inv(T4) * L * T4 ;
Mg = inv(T4) * M * T4 ;

Lg = M2T(Lg,'dot');
Mg = M2T(Mg,'dot');

%% GAUSSIAN INTEGRATION

% Retrieve Gauss points
[int_point,int_weight] = GaussPoints(2);

% Initialise K and B matrices
Ke  = zeros(ndsPerElem*dofN,ndsPerElem*dofN); 

for i = 1:numel(int_point)
     
    xi = int_point(i);
   
    for j = 1:numel(int_point)
        
        eta = int_point(j);
        
        % Shape function derivatives (Natural Coords):
        R = 1/4.*[(-1+eta) (1-eta) (1+eta) (-1-eta); 
                   (-1+xi) (-1-xi)  (1+xi)   (1-xi)];

        % Jacobian Matrix
        J =  R * elCoords;

        % Shape function derivatives (Global Coords):
        dN = inv(J)*R;

        % Displacement Gradient Relationship
        B = B_Matrix_2D(dN);

        % Displacement Gradient Vector
        GradU = B * q;

        % Assemble A1 matrix:
        A1 = [2 0 0 0;
              0 2 0 0;
              0 0 1 1;
              0 0 1 1];

        % Assemble A2 matrix:
        A2 = [GradU(1)   0        0      GradU(4);
                  0    GradU(2) GradU(3)   0     ;
                  0    GradU(4) GradU(1)   0     ;
              GradU(3)   0        0      GradU(2)];

        % Assemble A3 matrix:
        A3 = [2*GradU(1)     0          0      2*GradU(4);
                  0      2*GradU(2) 2*GradU(3)     0     ;
              1*GradU(3) 1*GradU(4) 1*GradU(1) 1*GradU(2);
              1*GradU(3) 1*GradU(4) 1*GradU(1) 1*GradU(2)];

        % Right Cauchy-Green deformation vector:
        Cvec = [1;1;0;0] + (A1 + A2)*GradU;
        
        % Right Cauchy-Green deformation tensor:
        C = V2M(Cvec);
        
        % Deformation gradient tensor:
        F = I + V2M(GradU);
        
        %% NEED TO CHANGE WHAT IS WRITTEN BELOW!
        % Simplified form of strain derivatives (small displacement)
        
        [V,E] = eig(C);
        normE = norm(diag(E));
        tol = 10^-3*normE;

        % Difference between eigenvalues
        d12 = abs(E(1,1)-E(2,2));

        % Offset one eigenvalue if two are equal (make them distinct)
        if d12<tol
            if E(1,1) >= E(2,2);  E(2,2) = E(2,2) - tol; end
            if E(2,2) >= E(1,1);  E(1,1) = E(1,1) - tol; end
        end
        Ct = V*E*inv(V);
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
        
        
        CC_mu = zeros(2,2,2,2); Sbar_M = zeros(2,2); SM = Sbar_M;
        for k = 1:num_param
            ni = nset(k);
            ci = cset(k);
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
        
        Cg = 2*(CC_lam + CC_mu + Bonus_);     % factor of 2 because CC = 2*(dS/dC);
        
        S = SL + SM;
        
        %% THIS HAS CHANGED BASED ON DERIVATION
        % Element tangent stiffness
        K = (B' * T2M(circledot22(I,S) + SinCon24(F,Transpose4(SinCon24(F,Cg),'dot'))) * B) .* det(J);
%         K = (B' * (A1+A3)' * Cg * (A1+A2) * B) .* det(J) ;
        
        Ke = Ke + (int_weight(i) * int_weight(j)) .* K ;

    end
    hello = 1;
    
    
end
