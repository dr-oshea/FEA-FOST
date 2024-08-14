function [DispN,DispO,DiffN,DiffO,Phi,Beta,IStage,PLambdaN,PhiN,...
    PhiO,NIter] = ProcessMod(DispN,DispO,DiffN,DiffO,Phi,Beta,NSF,IStage,...
    PLambdaN,PhiN,PhiO,NIter)

DiffO=DiffN;
DiffN = DispN - DispO;
Sum1 = DiffO'*DiffO;
Sum2 = DiffN'*DiffO;
PhiN=DiffN;
DispO=DispN;
PLambdaO = PLambdaN;
if (Sum1 == 0)
    % Do Nothing
else
    PLambdaN = Sum2 * sign(Sum2) / Sum1;
    LambdaTol = abs(PLambdaN - PLambdaO);

    PhiSign = sign(Sum2);
    SSS = (Sum1)^0.5;
    SSS = abs(SSS) * PhiSign;
    PH1 = DiffN(1) / SSS;
    if (PH1 < 0)
        PhiSign = -1;
    end
    PhiN = PhiN*PhiSign/SSS;
    PhiTol = 0.0;
    for i=1:NSF
        PhiTol = PhiTol + (PhiN(i) - PhiO(i))^2;
    end    
    PhiO=PhiN;
    if (PhiTol > 0.01 || LambdaTol > 0.01 || IStage >= 100 || abs(1 - PLambdaN) < 0.01) % sqrt(PhiTol)
        %Do Nothing
    else
        IStage = IStage+1;
        Beta(IStage) = PLambdaN / (1.0 - PLambdaN);
        for i=1:NSF
            Phi(i,IStage) = PhiN(i);
        end
    end
    fprintf('\n         *  Lambda = %.5f  ,  norm(Phi) = %.5f',PLambdaN,norm(PhiN));
    
end


