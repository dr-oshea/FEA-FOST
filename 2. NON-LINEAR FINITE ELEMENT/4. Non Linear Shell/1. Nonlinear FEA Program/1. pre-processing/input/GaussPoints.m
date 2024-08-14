function [P,W]=GaussPoints(num)
    
    if num == 1

        P = 0.000;
        W = 2.000;
        
    elseif num == 2
        
        P = [-0.57735026919  0.57735026919];
        W = [ 1.00000000000  1.00000000000];
        
    elseif num == 3
        
        P = [-0.7745966692  0.0000000000  0.7745966692];
        W = [ 0.5555555556  0.8888888889  0.5555555556];
        
    elseif num == 4
        
        P = [-0.8611363116 -0.3399810436  0.3399810436  0.8611363116];
        W = [ 0.3478548451  0.6521451549  0.6521451549  0.3478548451];
        
    elseif num == 5
        
        P = [-0.9061798459 -0.5384693101  0.0000000000  0.5384693101  0.9061798459];
        W = [ 0.2369268851  0.4786286705  0.5688888889  0.4786286705  0.2369268851];
        
    elseif num == 6
        
        P = [-0.9324695142 -0.6612093865 -0.2386191861  0.2386191861  0.6612093865  0.9324695142];
        W = [ 0.1713244924  0.3607615730  0.4679139346  0.4679139346  0.3607615730  0.1713244924]; 
        
    elseif num == 7
        
        P = [-0.9491079123 -0.7415311856 -0.4058451514  0.0000000000  0.4058451514  0.7415311856  0.9491079123];
        W = [ 0.1294849662  0.2797053915  0.3818300505  0.4179591837  0.3818300505  0.2797053915  0.1294849662]; 
        
    elseif num == 30
        
        W = [0.102852653000000,0.102852653000000,0.101762390000000,0.101762390000000,0.0995934210000000,0.0995934210000000,0.0963687370000000,0.0963687370000000,0.0921225220000000,0.0921225220000000,0.0868997870000000,0.0868997870000000,0.0807558950000000,0.0807558950000000,0.0737559750000000,0.0737559750000000,0.0659742300000000,0.0659742300000000,0.0574931560000000,0.0574931560000000,0.0484026730000000,0.0484026730000000,0.0387991930000000,0.0387991930000000,0.0287847080000000,0.0287847080000000,0.0184664680000000,0.0184664680000000,0.00796819200000000,0.00796819200000000];
        P = [-0.0514718430000000,0.0514718430000000,-0.153869914000000,0.153869914000000,-0.254636926000000,0.254636926000000,-0.352704726000000,0.352704726000000,-0.447033770000000,0.447033770000000,-0.536624148000000,0.536624148000000,-0.620526183000000,0.620526183000000,-0.697850495000000,0.697850495000000,-0.767777432000000,0.767777432000000,-0.829565762000000,0.829565762000000,-0.882560536000000,0.882560536000000,-0.926200047000000,0.926200047000000,-0.960021865000000,0.960021865000000,-0.983668123000000,0.983668123000000,-0.996893484000000,0.996893484000000];
        
    end
    
end