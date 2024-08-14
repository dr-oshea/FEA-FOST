function plotLoadShaded(input)

% -------------------------------------------------------------------------

% plot x,y,z displacements
% plot L,T,M displacements
% plot x,y,z 
% plot L,T,M stresses 

% I believe this conains all displacements at all nodes.
% resp.static.D = Q;
% dof_ = size(input.ND,1)*dofN;

% set up plots
NF = 2;

figure;
axis equal; axis off; hold on; view(3);
title('Displacements from static loads');
warning off

maxX = max(input.ND(:,2));   minX = min(input.ND(:,2));
maxY = max(input.ND(:,3));   minY = min(input.ND(:,3));
maxZ = max(input.ND(:,4));   minZ = min(input.ND(:,4));
labx = (maxX / 11); laby = (maxY / 11); labz = (maxZ / 11);
labx = max([labx laby labz]); laby = labx;  labz=labx;

% deN = max([max(abs(resp.static.D(1:dofN:dof_(1)))) ...
%       max(abs(resp.static.D(2:dofN:dof_(1)))) ...
%       max(abs(resp.static.D(3:dofN:dof_(1))))]);
% dx = labx *resp.static.D(1:dofN:dof_(1))./deN;
% dy = laby *resp.static.D(2:dofN:dof_(1))./deN;
% dz = labz *resp.static.D(3:dofN:dof_(1))./deN;
% 
ND_d = input.ND;
ND_d(:,2) = input.ND(:,2);
ND_d(:,3) = input.ND(:,3);
ND_d(:,4) = input.ND(:,4);
% 
[qu] = fem_3_surfBrick(input.ND(:,2:end),input.EL(:,3:10));
ndPlot = unique(qu);
% 
%plot_values = zeros(size(input.ND,1),1);

plot_values = linspace(10,10,size(input.ND,1))'; %making linspace numbers the same keeps the colours constant

% plot_values = resp.static.D(1:dofN:dof_(1));
% t = 'X-Displacements';

plot3(ND_d(ndPlot,2),ND_d(ndPlot,3),ND_d(ndPlot,4),'b.'); hold on;
for i=1:size(qu,1)
    patch([ND_d(qu(i,1),2) ND_d(qu(i,2),2) ND_d(qu(i,3),2) ...
           ND_d(qu(i,4),2)], ...
          [ND_d(qu(i,1),3) ND_d(qu(i,2),3) ND_d(qu(i,3),3) ...
           ND_d(qu(i,4),3)], ...
          [ND_d(qu(i,1),4) ND_d(qu(i,2),4) ND_d(qu(i,3),4) ...
           ND_d(qu(i,4),4)], ...
           [plot_values(qu(i,1)) plot_values(qu(i,2)) plot_values(qu(i,3)) ...
           plot_values(qu(i,4))]);
end

for i = 1:size(input.LOAD)
     quiver3(input.ND(input.LOAD(i,1),2),input.ND(input.LOAD(i,1),3),input.ND(input.LOAD(i,1),4),-input.LOAD(i,2),...
         -input.LOAD(i,3),-input.LOAD(i,4),0.0005,'ShowArrowHead','off',Color='r')
end

axis equal; axis off;

hold on; material dull; view(3); %colormap(hsv);
       camlight; lighting none; shading faceted;
% set(NF,'name',['  Deformed shape. MAX(x) = '  num2str(max(ND_d(:,2))) ...
%       '   MAX(y) = ' num2str(max(ND_d(:,3))) '   MAX(z) = ' ...
%       num2str(max(ND_d(:,4)))],'NumberTitle','off');
rotate3d(gca);set(gcf,'Pointer','arrow');
%colorbar
hold off
title('Discretisation Shaded')
end
