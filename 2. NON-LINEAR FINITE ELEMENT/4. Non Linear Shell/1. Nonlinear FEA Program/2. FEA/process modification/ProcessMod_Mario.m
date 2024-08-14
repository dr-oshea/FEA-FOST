function [QMN,QM,DiffN,DiffO,Phi,Beta,IStage,TLamdN,PhiN,...
    PhiO] = ProcessMod_Mario(QMN,QM,DiffN,DiffO,Phi,Beta,NSF,IStage,...
    TLamdN,PhiN,PhiO,alpha,maxMOD)

Sum1 = 0;
Sum2 = 0;
DiffO = DiffN;
DiffN = QMN - QM;
PhiN = DiffN;
QM = QMN;
for i = 1:NSF
    Sum1 = Sum1 + DiffO(i)*DiffO(i);
    Sum2 = Sum2 + DiffN(i)*DiffO(i);
end

TLamdO = TLamdN;

if Sum1 == 0; return; end

TLamdN = Sum2/Sum1;

ACC = abs(TLamdN-TLamdO);
TTL = abs(1.0 - TLamdN);

SSS = sqrt(Sum1);
SSS = SSS * sign(Sum2);

PSign = 1.0;

Ph1 = DiffN(1)/SSS;
if Ph1 < 0
    PSign = -1;
end

ACCPH = 0;
for i = 1:NSF
    PhiN(i) = PSign*PhiN(i)/SSS;
    ACCPH = ACCPH + (PhiN(i)-PhiO(i))^2;
end
PhiO = PhiN;

ACCPH = sqrt(ACCPH);

% Check tolerances
if ACC > alpha || TTL < 0.005 || ACCPH > alpha || IStage > maxMOD
    return; 
end
    
IStage = IStage + 1;
Beta(IStage) = TLamdN/(1.0-TLamdN);
for i = 1:NSF
    Phi(i,IStage) = PhiN(i);
end

end