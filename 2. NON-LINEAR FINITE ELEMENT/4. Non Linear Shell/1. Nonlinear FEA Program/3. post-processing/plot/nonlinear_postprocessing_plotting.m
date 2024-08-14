function nonlinear_postprocessing_plotting(input, Q, S, dofN) 

% -------------------------------------------------------------------------

% plot x,y,z displacements
% plot L,T,M displacements
% plot x,y,z 
% plot L,T,M stresses 

% I believe this conains all displacements at all nodes.
resp.static.D = Q;
dof_ = size(input.ND,1)*dofN;

% set up plots
num_plots = 10;


NF = 2;

for kk = 1:num_plots

figure;
axis equal; axis off; hold on; view(3);
title('Displacements from static loads');
warning off

maxX = max(input.ND(:,2));   minX = min(input.ND(:,2));
maxY = max(input.ND(:,3));   minY = min(input.ND(:,3));
maxZ = max(input.ND(:,4));   minZ = min(input.ND(:,4));
labx = (maxX / 11); laby = (maxY / 11); labz = (maxZ / 11);
labx = max([labx laby labz]); laby = labx;  labz=labx;

deN = max([max(abs(resp.static.D(1:dofN:dof_(1)))) ...
      max(abs(resp.static.D(2:dofN:dof_(1)))) ...
      max(abs(resp.static.D(3:dofN:dof_(1))))]);
dx = labx *resp.static.D(1:dofN:dof_(1))./deN;
dy = laby *resp.static.D(2:dofN:dof_(1))./deN;
dz = labz *resp.static.D(3:dofN:dof_(1))./deN;

ND_d = input.ND;
ND_d(:,2) = input.ND(:,2)+dx;
ND_d(:,3) = input.ND(:,3)+dy;
ND_d(:,4) = input.ND(:,4)+dz;

if 1
    [qu] = fem_3_surfBrick(input.ND(:,2:end),input.EL(:,3:10));
    ndPlot = unique(qu);
    
    if kk == 1
        plot_values = resp.static.D(1:dofN:dof_(1));
        t = 'X-Displacements';
    end

    if kk == 2
        plot_values = resp.static.D(2:dofN:dof_(1));
        t = 'Y-Displacements';
    end

    if kk == 3
        plot_values = resp.static.D(3:dofN:dof_(1));
        t = 'Z-Displacements';
    end

    if kk == 4
        plot_values = sqrt( (resp.static.D(1:dofN:dof_(1))).^2 + ...
       (resp.static.D(2:dofN:dof_(1))).^2 + (resp.static.D(3:dofN:dof_(1))).^2 );
        t = 'Total Resultant Displacement';
    end

    if kk == 5
        plot_values = S(:,1);
        t = 'FF-Stresses';
    end

    if kk == 6
        plot_values = S(:,2);
        t = 'SS-Stresses';
    end

    if kk == 7
        plot_values = S(:,3);
        t = 'NN-Stresses';
    end

    if kk == 8
        plot_values = S(:,4);
        t = 'SN-Stresses';
    end

    if kk == 9
        plot_values = S(:,5);
        t = 'FN-Stresses';
    end

    if kk == 10
        plot_values = S(:,6);
        t = 'FS-Stresses';
    end



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
    axis equal; axis off;
end

hold on; material dull; view(3); colormap(hsv);
       camlight; lighting none; % shading interp;
set(NF,'name',['  Deformed shape. MAX(x) = '  num2str(max(ND_d(:,2))) ...
      '   MAX(y) = ' num2str(max(ND_d(:,3))) '   MAX(z) = ' ...
      num2str(max(ND_d(:,4)))],'NumberTitle','off');
rotate3d(gca); set(gcf,'Pointer','arrow');
colorbar
hold off
title(t)

end

% SIGvonMises = zeros(length(SIG_main),1);
% for i=1:length(SIG_main)
%     SIGvonMises(i) = sqrt(0.5*((SIG_main(i,1)-SIG_main(i,2))^2 + ...
%         (SIG_main(i,2)-SIG_main(i,3))^2 + (SIG_main(i,1)-SIG_main(i,3)) ));
% end
% 
% figure(NF+1);
% axis equal; axis off; hold on; view(3);
% title('Stress field');
% KL = size(in_data.EL,1);
% com = 2;
% ts = 3;
% 
% if 1
%     for i=1:size(qu,1)
%         patch([ND_d(qu(i,1),2) ND_d(qu(i,2),2) ND_d(qu(i,3),2) ...
%                ND_d(qu(i,4),2)], ...
%               [ND_d(qu(i,1),3) ND_d(qu(i,2),3) ND_d(qu(i,3),3) ...
%                ND_d(qu(i,4),3)], ...
%               [ND_d(qu(i,1),4) ND_d(qu(i,2),4) ND_d(qu(i,3),4) ...
%                ND_d(qu(i,4),4)], ...
%               [SIGvonMises(qu(i,1)) SIGvonMises(qu(i,2)) ...
%                SIGvonMises(qu(i,3)) SIGvonMises(qu(i,4))]);
%     end
%     axis equal; axis off;
% end
% 
% hold on; material dull; view(3); colormap(hsv);
% [Smax,i] = max(SIG_main(:,1)); [Smin,k] = min(SIG_main(:,1));
% set(NF+1,'name',['  Max - node ' num2str(i) ': '  num2str(Smax) ...
%       '/' num2str(max(SIG_main(i,2))) '/' ...
%       num2str(max(SIG_main(i,3))) '. Min - node '...
%       num2str(k) ': ' num2str(Smin) '/' num2str(max(SIG_main(k,2))) '/' ...
%       num2str(max(SIG_main(k,3)))],'NumberTitle','off');
% colorbar('vert');
% 
% 
% maxX = max(in_data.ND(:,2));   minX = min(in_data.ND(:,2));
% maxY = max(in_data.ND(:,3));   minY = min(in_data.ND(:,3));
% maxZ = max(in_data.ND(:,4));   minZ = min(in_data.ND(:,4));
% labx = (maxX / 11); laby = (maxY / 11); labz = (maxZ / 11);
% labx = min([labx laby labz]); laby = labx;  labz=labx;
% 
% for i=1:size(in_data.CON)
%    node_i = find(in_data.ND(:,1)==in_data.CON(i,1));
%    if in_data.CON(i,2)==0
%       plot3([in_data.ND(node_i,2) in_data.ND(node_i,2)-labx], ...
%           [in_data.ND(node_i,3) in_data.ND(node_i,3)],...
%           [in_data.ND(node_i,4) in_data.ND(node_i,4)],'k-','LineWidth',2);
%           hold on;
%    end;
%    if in_data.CON(i,3)==0
%       plot3([in_data.ND(node_i,2) in_data.ND(node_i,2)], ...
%           [in_data.ND(node_i,3) in_data.ND(node_i,3)-laby],...
%           [in_data.ND(node_i,4) in_data.ND(node_i,4)],'k-','LineWidth',2);
%           hold on;
%    end;
%    if in_data.CON(i,4)==0
%       plot3([in_data.ND(node_i,2) in_data.ND(node_i,2)], ...
%           [in_data.ND(node_i,3) in_data.ND(node_i,3)], ...
%           [in_data.ND(node_i,4) in_data.ND(node_i,4)-labz], ...
%           'k-','LineWidth',2); hold on;
%    end;
%    if in_data.CON(i,4)==0 | in_data.CON(i,5)==0 | in_data.CON(i,6)==0
%       plot3([in_data.ND(node_i,2) in_data.ND(node_i,2)], ...
%           [in_data.ND(node_i,3) in_data.ND(node_i,3)],...
%           [in_data.ND(node_i,4) in_data.ND(node_i,4)],'ks','LineWidth',2);
%           hold on;
%    end;
% end;
% 
% rotate3d(gca); set(gcf,'Pointer','arrow');
% warning on
