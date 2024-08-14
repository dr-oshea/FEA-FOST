%% Improved Guess
function [X_new] = ImprovedGuess(Phi,Beta,XN,XO,mod,dofT)

% This script calculates an improves guess based on the process
% modification technique in Lawther (1980). It uses the eigenvalue of the
% iterative process to improve convergence characteristics

% Author:   Daniel O'Shea
% Created:  05 April 2018

% INPUTS:
% Phi = array storing Phi vector for each modified process
% Beta = array storing Beta scalar for each modified process
% XN = vector containing initial estimate of displacement
% XO = vector containing [?????]
% mod = total stages of modifications
% dofT = total degrees of freedom of model

% OUTPUTS:
% X_new = improved estimated displacement vector


%% ------------------------------------------------------------------------

% % For each stage of process modification
% for i = 1 : mod
%     Sum = 0;    % reset sum
%     
%     % Calculate first part of improved guess
%     for j = 1 : dofT    
%         Sum = Sum + Phi(j,i) * (XN(j) - XO(j));
%     end
%     
%     % Multiply by current beta
%     Sum = Sum * Beta(i);
%     
%     % Complete the improved guess
%     for j = 1 : dofT
%         XN(j) = XN(j) + Sum * Phi(j,i);
%     end 
%     
% end

% Vector form:
for i = 1 : mod
    
    Sum = Beta(i) * Phi(:,i)' * (XN - XO);
    
    % Complete the improved guess
    for j = 1 : dofT
        XN(j) = XN(j) + Sum * Phi(j,i);
    end 
    
end

X_new = XN;

