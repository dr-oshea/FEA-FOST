function [C0,L0,M0] = input_material(MatProps,PLANAR,type)

%%  input_material

% This script generates three orthotropic material tensors depending on the
% prescribed material properties from the input file. Works for 2D and 3D
% problems

% Author:   Daniel O'Shea
% Created:  16 March 2018

% INPUTS:
% MatProps = material properties in terms of E, v, G values

% OUTPUTS:
% C0 = classical fourth order material tensor flattenned to matrix
% L0 = classical fourth order Lambda material tensor
% M0 = classical fourth order Mu material tensor

% UPDATE LOG:
% 12/04/2018 - included support for a 3D element
% 13/04/2018 - included support for plane strain problems


%% ---------------------------------------------------------------------------

E11 = MatProps(1);
E22 = MatProps(2);
E33 = MatProps(3);
v23 = MatProps(4);
v13 = MatProps(5);
v12 = MatProps(6);
G23 = MatProps(7);
G13 = MatProps(8);
G12 = MatProps(9);

v21 = E22/E11*v12;
v31 = E33/E11*v13;
v32 = E33/E22*v23;

% Components of compliance matrix
S1111 = 1/E11;     S2222 = 1/E22;     S3333 = 1/E33;
S1122 = -v21/E22;  S3311 = -v13/E11;  S2233 = -v32/E33;
S2211 = -v12/E11;  S1133 = -v31/E33;  S3322 = -v23/E22;
S2323 = 1/(2*G23); S3131 = 1/(2*G13); S1212 = 1/(2*G12);
S3232 = 1/(2*G23); S1313 = 1/(2*G13); S2121 = 1/(2*G12);
    
% Compliance Matrix
S = [S1111 S1122 S1133   0     0     0     0     0     0   ;
     S2211 S2222 S2233   0     0     0     0     0     0   ;
     S3311 S3322 S3333   0     0     0     0     0     0   ;
       0     0     0   S2323   0     0     0     0     0   ;
       0     0     0     0   S3131   0     0     0     0   ;
       0     0     0     0     0   S1212   0     0     0   ;
       0     0     0     0     0     0   S3232   0     0   ;
       0     0     0     0     0     0     0   S1313   0   ;
       0     0     0     0     0     0     0     0   S2121];
   
% Stiffness Matrix
C = inv(S);

% Components of stiffness matrix
C1111 = C(1,1); C1122 = C(1,2); C1133 = C(1,3);
C2211 = C(2,1); C2222 = C(2,2); C2233 = C(2,3);
C3311 = C(3,1); C3322 = C(3,2); C3333 = C(3,3);
C2323 = C(4,4); C3131 = C(5,5); C1212 = C(6,6);
C3232 = C(7,7); C1313 = C(8,8); C2121 = C(9,9);

%% 2D Constitutive Law
if type == 1
    
    % Plane Stress
    if PLANAR == 1
        % Adjusted constants
        c1111 = C1111-C1133*C3311/C3333;
        c1122 = C1122-C1133*C2233/C3333;
        c2222 = C2222-C2233*C3322/C3333;
        c1212 = C1212;
        c2121 = C2121;
    
    % Plane Strain
    elseif PLANAR == 2
        % Adjusted constants
        s1111 = S1111-S1133*S3311/S3333;
        s1122 = S1122-S1133*S3322/S3333;
        s2211 = S2211-S3311*S2233/S3333;
        s2222 = S2222-S2233*S3322/S3333;
        s1212 = S1212;
        s2121 = S2121;
        
        S_ = [s1111 s1122   0     0   ;
              s2211 s2222   0     0   ;
                0     0   s1212   0   ;
                0     0     0   s2121];
         
        C_ = inv(S_);
        
        c1111 = C_(1,1); c2222 = C_(2,2); c1122 = C_(1,2); 
        c1212 = C_(3,3); c2121 = C_(4,4);
        
    end
        
    L1111 = c1122;
    L2222 = c1122;
    L1122 = c1122;

    M1111 = c1111-L1111;
    M2222 = c2222-L2222;
    M1212 = c1212;
    M2121 = c2121;

    L0 = [L1111 L1122 0 0;
          L1122 L2222 0 0;
            0     0   0 0;
            0     0   0 0];

    M0 = diag([M1111 M2222 M1212 M2121]);

end

%% 3D Constitutive Law
if type == 2
      
    L1111 = + C1122 + C1133 - C2233;
    L2222 = + C1122 - C1133 + C2233;
    L3333 = - C1122 + C1133 + C2233;
    L1122 = (1/2)*(L1111 + L2222);
    L1133 = (1/2)*(L1111 + L3333);
    L2233 = (1/2)*(L2222 + L3333);
    
    M1111 = C1111-L1111;
    M2222 = C2222-L2222;
    M3333 = C3333-L3333;
    M2323 = C2323;
    M1313 = C1313;
    M1212 = C1212;
    
    L0 = [L1111 L1122 L1133 0  0  0  0  0  0;
          L1122 L2222 L2233 0  0  0  0  0  0;
          L1133 L2233 L3333 0  0  0  0  0  0;
          0     0     0     0 0 0 0 0 0;
          0     0     0     0 0 0 0 0 0;
          0     0     0     0 0 0 0 0 0;
          0     0     0     0 0 0 0 0 0;
          0     0     0     0 0 0 0 0 0;
          0     0     0     0 0 0 0 0 0];
      
    M0 = diag([M1111; M2222; M3333; M2323; M1313; M1212; M2323; M1313; M1212]);
    
end


C0 = L0 + M0;


end

