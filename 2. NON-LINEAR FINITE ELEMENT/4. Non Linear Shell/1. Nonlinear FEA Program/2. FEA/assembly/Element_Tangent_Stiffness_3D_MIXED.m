function [KT, tvec] = Element_Tangent_Stiffness_3D_MIXED (q, elCoords_X, theta, model_data, ~)

% This script generates the element tangent stiffness matrix at a 
% given displacement field for a 3D brick element using a three field mixed
% finite deformation formulation. The integral is determined using Gaussian
% integration. The formulation here is from Zienciewicz and Taylor
% "The Finite Element Method for Solid and Structural Mechanics"

% Author:   Daniel O'Shea
% Created:  22 June 2023

% INPUTS:
% q = local vector of known displacements
% elCoords = nodal cordinates for local nodes
% theta = material orientation

% OUTPUTS:
% K0e = element stiffness matrix

%% ---------------------------------------------------------------------------
ndsPerElem = 8;
dofN = 3;
phi = deg2rad(theta);
I = Iden2(3);
Ivec = M2V(I);
IoI = circledot22(I,I);
IxI = circlecross22(I,I);

% Global Element Constitutive Law
% [~,~,T9,~,~] = Transform2(phi);

%% Updated coordinates:

Qmat = reshape(q',[3,8]).';
elCoords_x = elCoords_X + Qmat;

%% GAUSSIAN INTEGRATION

% Create means of returning C tensor and S tensor for various motions (have
% deformation gradient as an argument?)

% Calculate intermediate solutions for theta and p (integration)
%   K_tp, Pt, Pp
% Use integration to get the various stiffness matrices
%   K_uu, K_up, K_ut, K_tt
% Calculate KT

% Retrieve Gauss points
[int_point_1,int_weight_1] = GaussPoints(1);
[int_point_2,int_weight_2] = GaussPoints(2);

% Initialise K matrices
Pp = 0;
Pt = 0;
Kpt = 0; 

% Perform numeric integration for Kpt, Pp

for i = 1:numel(int_point_1)
     
    xi = int_point_1(i);
   
    for j = 1:numel(int_point_1)
        
        eta = int_point_1(j);
        
        for k = 1:numel(int_point_1)
        
            mu = int_point_1(k);
            
            % Shape functions for p and theta
            Np = 1;
            Nt = 1;

            % Shape function derivatives (Natural Coords) - BMATRIX:
            R = 1/8.*[-(1-eta)*(1-mu) (1-eta)*(1-mu) (1+eta)*(1-mu) -(1+eta)*(1-mu) -(1-eta)*(1+mu) (1-eta)*(1+mu) (1+eta)*(1+mu) -(1+eta)*(1+mu);...
               -(1-xi)*(1-mu) -(1+xi)*(1-mu) (1+xi)*(1-mu) (1-xi)*(1-mu) -(1-xi)*(1+mu) -(1+xi)*(1+mu) (1+xi)*(1+mu) (1-xi)*(1+mu);...
               -(1-xi)*(1-eta) -(1+xi)*(1-eta) -(1+xi)*(1+eta) -(1-xi)*(1+eta) (1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

            JX = R * elCoords_X;
            dNX = inv(JX)*R;
            BX = B_Matrix_3D(dNX);

            % Displacement Gradient Tensor
            GradU = V2M(BX * q);

            % Deformation gradient vector
            F = GradU + I;
            J = det(F);
            
            % Calculate Kpt, Pp
            weight = int_weight_1(i) * int_weight_1(j) * int_weight_1(k);
            Kpt = Kpt + weight .* (Nt.'*Np * det(JX)) ;
            Pp = Pp + weight .* (Np.'* J .* det(JX)) ;
        
        end
    end
end

% Partial solution of theta for element:
tvec = Kpt\Pp;

% Perform numeric integration for Kpt, Pp

for i = 1:numel(int_point_1)
     
    xi = int_point_1(i);
   
    for j = 1:numel(int_point_1)
        
        eta = int_point_1(j);
        
        for k = 1:numel(int_point_1)
        
            mu = int_point_1(k);
            
            % Shape functions for theta
            Nt = 1;

            theta = Nt*tvec;

            % Shape function derivatives (Natural Coords) - BMATRIX:
            R = 1/8.*[-(1-eta)*(1-mu) (1-eta)*(1-mu) (1+eta)*(1-mu) -(1+eta)*(1-mu) -(1-eta)*(1+mu) (1-eta)*(1+mu) (1+eta)*(1+mu) -(1+eta)*(1+mu);...
               -(1-xi)*(1-mu) -(1+xi)*(1-mu) (1+xi)*(1-mu) (1-xi)*(1-mu) -(1-xi)*(1+mu) -(1+xi)*(1+mu) (1+xi)*(1+mu) (1-xi)*(1+mu);...
               -(1-xi)*(1-eta) -(1+xi)*(1-eta) -(1+xi)*(1+eta) -(1-xi)*(1+eta) (1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

            JX = R * elCoords_X;
            dNX = inv(JX) * R;
            BX = B_Matrix_3D(dNX);

            % Displacement Gradient Tensor
            GradU = V2M(BX * q);

            % Deformation gradient vector
            F = GradU + I;
            J = det(F);
            F_mod = (theta / J)^(1/3)*F;
            C_mod = F_mod.'*F_mod;

            % Get PK2 stress tensor, and tangent stiffness tensor           
            [S_mod,~] = Get_S_CC(F_mod, model_data);
            
            % Get Cauchy stress tensor and current config elasticity tensor
            Pbar = 1/3*theta^(-1)*(DblCon22(S_mod,C_mod));
            
            % Calculate Kpt, Pp
            weight = int_weight_1(i) * int_weight_1(j) * int_weight_1(k);
            Pt = Pt + weight .* (Nt.'* Pbar .* det(JX)) ;
        
        end
    end
end

% Partial solution of p for element:
pvec = Kpt\Pt;

% debugging
if abs(tvec) > 1.001
    hello = 1;
end


%% Numerical integration to Determine tangent stiffness

Kuu = zeros(ndsPerElem*dofN,ndsPerElem*dofN); 
Kup = zeros(ndsPerElem*dofN,1); 
Kut = zeros(ndsPerElem*dofN,1); 
Ktt = 0; 

Kg_ = 0;

% Assemble matrix for fourth-order symmetry projection tensor:
S4 = (1/2)*(circledot22(I,I) + circlestar22(I,I));


for i = 1:numel(int_point_2)
     
    xi = int_point_2(i);
   
    for j = 1:numel(int_point_2)
        
        eta = int_point_2(j);
        
        for k = 1:numel(int_point_2)
        
            mu = int_point_2(k);
            
            % Shape functions for p and theta
            Np = 1;
            Nt = 1;
            theta = Nt*tvec;
            p = Np*pvec;

            % Shape function derivatives (Natural Coords) - BMATRIX:
            R = 1/8.*[-(1-eta)*(1-mu) (1-eta)*(1-mu) (1+eta)*(1-mu) -(1+eta)*(1-mu) -(1-eta)*(1+mu) (1-eta)*(1+mu) (1+eta)*(1+mu) -(1+eta)*(1+mu);...
               -(1-xi)*(1-mu) -(1+xi)*(1-mu) (1+xi)*(1-mu) (1-xi)*(1-mu) -(1-xi)*(1+mu) -(1+xi)*(1+mu) (1+xi)*(1+mu) (1-xi)*(1+mu);...
               -(1-xi)*(1-eta) -(1+xi)*(1-eta) -(1+xi)*(1+eta) -(1-xi)*(1+eta) (1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

            % Jx =  R * elCoords_x;
            % dNx = inv(Jx) * R;
            % Bx = Bs_Matrix_3D(dNx);

            JX = R * elCoords_X;
            dNX = inv(JX) * R;
            BX = B_Matrix_3D(dNX);

            % Displacement Gradient Tensor
            GradU = V2M(BX * q);

            % Deformation gradient vector
            F = GradU + I;
            C = F.'*F;
            J = det(F);
            F_mod = (theta / J)^(1/3)*F;
            C_mod = F_mod.'*F_mod;
            invC = inv(C);

            % Get PK2 stress tensor, and tangent stiffness tensor           
            [S_mod,CC_mod] = Get_S_CC(F, model_data);
            Pbar = 1/3*theta^(-1)*DblCon22(S_mod,C_mod);

            % Get projection tensors
            PD = IoI - 1/3*(circlecross22(invC,C));
            PDt = IoI - 1/3*(circlecross22(C,invC));
            Ps = circledot22(invC,invC) - (1/3)*circlecross22(invC,invC);
            S4 = (1/2)*(circledot22(I,I) + circlestar22(I,I));
            
            % Assemble (I odot GradU) tensor, convert to matrix
            A1 = DblCon44(S4,circledot22(F.',I));
            A1t = Transpose4(A1,'dot');

            % Get 'complete' PK2 and elasticity tensors
            S_iso = (theta/J)^(2/3)*DblCon42(PD,S_mod);
            S_vol = p*J*invC;
            S_full = S_iso + S_vol;

            CC_iso = (theta/J)^(4/3)*DblCon44(PD,DblCon44(CC_mod,PDt)) + 2*theta*Pbar*Ps - (2/3)*(theta/J)^(4/3)*(circlecross22(DblCon42(PD,S_mod),inv(C)) + circlecross22(inv(C),DblCon24(S_mod,PDt)));
            CC_vol = (p*J)*circlecross22(inv(C),inv(C))-2*(p*J)*circledot22(inv(C),inv(C));
            CC_full = CC_vol + CC_iso;

            Duu = circledot22(I,S_full) + DblCon44(A1t,DblCon44(CC_full,A1));
            Dut = 1/3*theta^(-1)*(theta/J)^(2/3)*DblCon42(A1,DblCon42(PD,DblCon42(CC_mod,C_mod)+2*DblCon42(PD,S_mod)));
            Dup = J*DblCon42(A1t,inv(C));
            Dtt = 1/3*theta^(-1)*(1/3*DblCon22(C_mod,DblCon42(CC_mod,C_mod))-Pbar);
            

            % duu = DblCon44(Pd, DblCon44(cc_mod, Pd)) - 2/3*(circlecross22(sigma_dev,I) + circlecross22(I, sigma_dev)) - 2*(pmod - pbar)*IoI + (pmod - 2/3*pbar)*IxI;
            % dut = 2/3*sigma_dev + 1/3*DblCon42(Pd,DblCon42(cc_mod, I));
            % dtt = 1/theta*(1/9*DblCon22(I,DblCon42(cc_mod, I)) - 1/3*pbar);
            % kappa = model_data(2);

            % [~,Cvol] = Get_Theta_Deriv(model_data, t);

            % Calcualte element tangent stiffness in global coords
            Kuu_ = (BX' * T2M(Duu) * BX) .* det(JX);
            Kup_ = (BX' * M2V(Dup) * Np) .* det(JX);
            Kut_ = (BX' * M2V(Dut) * Nt) .* det(JX);
            Ktt_ = (Nt' * Nt .* Dtt).* det(JX);

            % Calculate geometric stiffness
            Kg_ = (BX' * T2M(circledot22(I, S_mod)) * BX) .* det(JX);
            
            weight = int_weight_2(i) * int_weight_2(j) * int_weight_2(k);
            Kuu = Kuu + weight .* (Kuu_ + Kg_);
            Kup = Kup + weight .* Kup_ ;
            Kut = Kut + weight .* Kut_ ;
            Ktt = Ktt + weight .* Ktt_ ;
        
        end
    end
    hello = 1;
end

Ktp = Kpt;
Kpu = Kup.';
Ktu = Kut.';

KT = Kuu + Kut*inv(Kpt)*Kpu + Kup*(Kpt\Ktu) + Kup*inv(Ktp)*Ktt*inv(Kpt)*Kpu;
