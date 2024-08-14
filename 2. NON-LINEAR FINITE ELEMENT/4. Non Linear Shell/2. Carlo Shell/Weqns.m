%% DIFFERENTIATION OF STRAIN ENERGY WITH RESPECT TO U
clear
clc
disp('    Beginning script: Determining Equations from Strain Energy')

% Define independent parameters of strain energy
syms u11 u12 u13 u22 u23 u33
assume(u11, 'real'); assume(u12, 'real'); assume(u13, 'real');
assume(u22, 'real'); assume(u23, 'real'); assume(u33, 'real');
syms L1111 L2222 L3333
syms M1111 M2222 M3333 M2323 M1313 M1212

Uvec = [u11; u22; u33; u23; u13; u12];
Lvec = [L1111; L2222; L3333];
Mvec = [M1111; M2222; M3333; M2323; M1313; M1212];
par = [Uvec; Lvec; Mvec];

Z = symfun(0,par);

disp('    Unknown parameters defined.')

% Establish random variables for testing purposes:
u11n = rand; u12n = rand; u13n = rand;
u22n = rand; u23n = rand; u33n = rand;

L1111n = rand; L2222n = rand; L3333n = rand; 
M1111n = rand; M2222n = rand; M3333n = rand;
M2323n = rand; M1313n = rand; M1212n = rand;

%% Stretch Tensor

% Symmetric Right Stretch Tensor U
U = [u11 u12 u13;
     u12 u22 u23;
     u13 u23 u33];
 
disp('    Right stretch tensor U defined.')

%% Eigendecomposition
disp('    Commencing Eigendecomposition of U.')

[Vec,Val] = eig(U);

% Define eigenvalues as functions of components of U
eU11 = symfun(Val(1,1),par);
eU22 = symfun(Val(2,2),par);
eU33 = symfun(Val(3,3),par);

% Diagonal Matrix of Eigenvalue, defined as function of components of U
eU = cell(3,3);
eU{1}{1} = real(eU11);
eU{2}{2} = real(eU22);
eU{3}{3} = real(eU33);
eU{1}{2} = Z; eU{2}{1} = Z;
eU{1}{3} = Z; eU{3}{1} = Z;
eU{2}{3} = Z; eU{3}{2} = Z;

% Eigenvector components:
Q11 = symfun(Vec(1,1),par); Q12 = symfun(Vec(1,2),par); Q13 = symfun(Vec(1,3),par);
Q21 = symfun(Vec(2,1),par); Q22 = symfun(Vec(2,2),par); Q23 = symfun(Vec(2,3),par);
Q31 = symfun(Vec(3,1),par); Q32 = symfun(Vec(3,2),par); Q33 = symfun(Vec(3,3),par);

Q_1 = {Q11; Q12; Q13};
Q_2 = {Q21; Q22; Q23};
Q_3 = {Q31; Q23; Q33};

Q = {Q_1 Q_2 Q_3};
Qt = transpose(Q);

% Notation used to store functions in Q:
% Q{i}{j} == Qij
disp('    Eigendecomposition of U complete.')

%% Fourth order eigenvector tensors:
disp('    Establishing fourth order eigenvector tensors...')

n = 3; % dimensions
QQ = cell(n,n,n,n);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                QQ{i}{l}{j}{k} = Q{i}{j}*Qt{k}{l};
            end
        end
    end
end


QQt = cell(n,n,n,n);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                QQt{i}{l}{j}{k} = Qt{i}{j}*Q{k}{l};
            end
        end
    end
end

disp('    Fourth order Q tensors generated.')

%% Material tensors

L1122 = (1/2)*(L1111+L2222);
L1133 = (1/2)*(L1111+L3333);
L2233 = (1/2)*(L2222+L3333);

L11 = symfun(L1111,par);
L22 = symfun(L2222,par);
L33 = symfun(L3333,par);
L12 = symfun(L1122,par);
L13 = symfun(L1133,par);
L23 = symfun(L2233,par);
Z = symfun(0,par);

% Matrix to Tensor
voigt = [1 6 8;
         9 2 4;
         5 7 3];
Lambda = cell(3,3,3,3);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                Lambda{i}{j}{k}{l}=Z;
            end
        end
    end
end      
Lambda{1}{1}{1}{1}=L11;
Lambda{2}{2}{2}{2}=L22;
Lambda{3}{3}{3}{3}=L33;
Lambda{2}{2}{3}{3}=L23;
Lambda{3}{3}{1}{1}=L13;
Lambda{1}{1}{2}{2}=L12;
Lambda{2}{2}{3}{3}=L23;
Lambda{1}{1}{3}{3}=L13;
Lambda{2}{2}{1}{1}=L12;
disp('    Lambda material tensor defined.')

% Random numerical tensor (for testing purposes):
L_num = cell(3,3,3,3);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                L_num{i}{j}{k}{l}=double(Lambda{i}{j}{k}{l}(u11n, u12n, u13n, u22n, u23n, u33n, L1111n, L2222n, L3333n, M1111n, M2222n, M3333n, M2323n, M1313n, M1212n));
            end
        end
    end
end 

%% Transformed Lambda Tensor
disp('    Commencing Transformation of Lambda...')

% Lhat = QQ : B
B = cell(3,3,3,3);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                B{i}{j}{k}{l}=Z;
            end
        end
    end
end
disp('    Initialised B cell.')

for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                for m=1:n
                    for n=1:n
                        B{i}{j}{k}{l} = B{i}{j}{k}{l} + Lambda{i}{j}{m}{n} * QQ{m}{n}{k}{l};
                    end
                end
            end
        end
    end
end
disp('    Computed Lambda * Q4.')

Lhat = cell(3,3,3,3);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                Lhat{i}{j}{k}{l}=Z;
            end
        end
    end
end
disp('    Initialising Transformed Lambda Cell Structure.')

for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                for m=1:n
                    for n=1:n
                        Lhat{i}{j}{k}{l} = Lhat{i}{j}{k}{l} + QQt{i}{j}{m}{n} * B{m}{n}{k}{l};
                    end
                end
            end
        end
    end
end
disp('    Computed Transformed Lambda.')

% Random Numeric testing:
Lh_num = cell(3,3,3,3);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                Lh_num{i}{j}{k}{l}=double(Lhat{i}{j}{k}{l}(u11n, u12n, u13n, u22n, u23n, u33n, L1111n, L2222n, L3333n, M1111n, M2222n, M3333n, M2323n, M1313n, M1212n));
            end
        end
    end
end 


%% Strain Energy

WL_h = 0;
for i = 1:n
    for j = 1:n
        for k = 1:n
            for l = 1:n
%                 WL = WL + (1/2)*Lambda{i}{i}{j}{j}*log(U{i}{i})*log(U{j}{j});
%                 WL_n = WL_n + (1/2)*L_num{i}{i}{j}{j}*log(U_num(i,i))*log(U_num(j,j));
                WL_h = WL_h + (1/2)*Lhat{i}{i}{j}{j}*log(eU{i}{i})*log(eU{j}{j});
%                 WL_h_n = WL_h_n + (1/2)*Lh_num{i}{i}{j}{j}*log(eU_num{i}{i})*log(eU_num{j}{j});
            end
        end
    end
end






