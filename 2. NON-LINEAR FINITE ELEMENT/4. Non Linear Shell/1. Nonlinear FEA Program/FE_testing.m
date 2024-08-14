%% Test Uniaxial Tension
div = 50;

% Material Properties
E = 1.43;
v = 0.49;
G = E/(2*(1+v));

E11 = E;
E22 = E;
E33 = E;
v23 = v;
v13 = v;
v12 = v;
G23 = G;
G13 = G;
G12 = G;
MatProps = [E11 E22 E33 v23 v13 v12 G23 G13 G12];
theta = 0;    % degrees

n1 = 0.75;
n2 = -0.703;
n3 = 2.628;
% nset = [n1];
nset = [n1; n2; n3];

% Mu additional constants
c2 = 0.091;
c3 = 0.0004;
cset = [1-c2-c3; c2; c3];

lamO = 1;
lamF = 1.2;

lamvec = lamO : (lamF-lamO)/div : lamF;

for i = 1 : numel(lamvec)
    lambda1 = lamvec(i);
    phi = deg2rad(theta);
    if lambda1 == 5.3;
        hello = 1;
    end
    
    if type == 1
        [S,F,W] = Uniaxial_StressTensor_2D(lambda1,MatProps,phi,nset,cset);
    elseif type == 2
        [S,F,W] = Uniaxial_StressTensor(lambda1,MatProps,phi,nset,cset);
    end
    
    tau = F * S * F';
    P = F * S';
    
    yplot(i) = tau(1,1);
    xplot(i) = lambda1;
    xplot2(i) = F(2,2);
    pplot(i) = P(1,1);
    
end

figure
hold on
plot(xplot,yplot,'-x','DisplayName','Hyperelastic Model')
plot(lamx,tauxx,'-o','DisplayName','FE Model');
ylim([0 max(max(yplot),max(tauxx))])
xlim([lamO max(max(xplot),max(lamx))])
legend('show')

figure
hold on
plot(xplot,pplot,'-x','DisplayName','Hyperelastic Model')
plot(lamx,Pxx,'-o','DisplayName','FE Model');
ylim([0 max(max(pplot),max(Pxx))])
xlim([lamO max(max(xplot),max(lamx))])
legend('show')

figure
hold on
plot(xplot,xplot2,'-x','DisplayName','Hyperelastic Model')
plot(lamx,lamy,'-o','DisplayName','FE Model');
ylim([min(min(xplot2),min(lamy)) 1])
xlim([lamO max(max(xplot),max(lamx))])
legend('show')