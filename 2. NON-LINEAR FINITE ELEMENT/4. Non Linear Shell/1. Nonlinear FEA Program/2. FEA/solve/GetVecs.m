%% Function to determine vectors
function [V,W] = GetVecs(i,DelQ,R,L)
    
% i = i+1;

DQ1_ = DelQ(:,i-1); DQ1 = DQ1_(L);
DQi_ = DelQ(:,i);   DQi = DQi_(L);

R1_ = R(:,i-1); R1 = R1_(L);
Ri_ = R(:,i);   Ri = Ri_(L);

gamma = R1-Ri;

% Determining V
numer = DQ1.'*gamma;
denom = DQi.'*R1;

bracket = 1-numer/denom;

V = bracket*R1 - Ri;

% Determining W
W = 1/numer * DQ1;

    
end