function plotFibre(input,PND,type)

% Adapted from Antony Ziecenco FEM Toolbox. Plots the problem
% discretisation and boundary conditions

%% ------------------------------------------------------------------------

NF = 110;

hold off;
%% Plot 3D brick
plot3(input.ND(:,2),input.ND(:,3),input.ND(:,4),'y.');
axis equal; axis off; view(3); hold on;

for i=1:size(input.EL)
    node1 = find(input.ND(:,1)==input.EL(i,3));
    node2 = find(input.ND(:,1)==input.EL(i,4));
    node3 = find(input.ND(:,1)==input.EL(i,5));
    node4 = find(input.ND(:,1)==input.EL(i,6));
    node5 = find(input.ND(:,1)==input.EL(i,7));
    node6 = find(input.ND(:,1)==input.EL(i,8));
    node7 = find(input.ND(:,1)==input.EL(i,9));
    node8 = find(input.ND(:,1)==input.EL(i,10));

    plot3([input.ND(node1,2) input.ND(node2,2) input.ND(node3,2) input.ND(node4,2) input.ND(node1,2) ...
          input.ND(node5,2) input.ND(node6,2) input.ND(node7,2) input.ND(node8,2) input.ND(node5,2)],...
       [input.ND(node1,3) input.ND(node2,3) input.ND(node3,3) input.ND(node4,3) input.ND(node1,3) ...
          input.ND(node5,3) input.ND(node6,3) input.ND(node7,3) input.ND(node8,3) input.ND(node5,3)],...
       [input.ND(node1,4) input.ND(node2,4) input.ND(node3,4) input.ND(node4,4) input.ND(node1,4) ...
          input.ND(node5,4) input.ND(node6,4) input.ND(node7,4) input.ND(node8,4) input.ND(node5,4)],...
          'Color',[0.4 0.1 0.7],'LineWidth',1);    

    plot3([input.ND(node2,2) input.ND(node6,2) input.ND(node7,2) input.ND(node3,2) input.ND(node2,2) ...
          input.ND(node1,2) input.ND(node5,2) input.ND(node8,2) input.ND(node4,2) input.ND(node1,2)],...
       [input.ND(node2,3) input.ND(node6,3) input.ND(node7,3) input.ND(node3,3) input.ND(node2,3) ...
          input.ND(node1,3) input.ND(node5,3) input.ND(node8,3) input.ND(node4,3) input.ND(node1,3)],...
       [input.ND(node2,4) input.ND(node6,4) input.ND(node7,4) input.ND(node3,4) input.ND(node2,4) ...
         input.ND(node1,4) input.ND(node5,4) input.ND(node8,4) input.ND(node4,4) input.ND(node1,4)],...
         'Color',[0.4 0.1 0.7],'LineWidth',1);
%plotting fibre vectors    
    element_centre_x_values = [input.ND(node1,2),input.ND(node2,2),input.ND(node3,2),input.ND(node4,2),...
        input.ND(node5,2),input.ND(node6,2),input.ND(node7,2),input.ND(node8,2)];
    element_centre_y_values = [input.ND(node1,3),input.ND(node2,3),input.ND(node3,3),input.ND(node4,3),...
        input.ND(node5,3),input.ND(node6,3),input.ND(node7,3),input.ND(node8,3)];
    element_centre_z_values = [input.ND(node1,4),input.ND(node2,4),input.ND(node3,4),input.ND(node4,4),...
        input.ND(node5,4),input.ND(node6,4),input.ND(node7,4),input.ND(node8,4)];
    element_centre_x = mean(element_centre_x_values);
    element_centre_y = mean(element_centre_y_values);
    element_centre_z = mean(element_centre_z_values);
    fibre_U = input.FIBRES.f(i,1);
    fibre_V = input.FIBRES.f(i,2);
    fibre_W = input.FIBRES.f(i,3);
    fibre_colour = rand(1,3);
    quiver3(element_centre_x,element_centre_y,element_centre_z,fibre_U,fibre_V,fibre_W,5,'ShowArrowHead','off',Color= fibre_colour, LineWidth=2);
    %quiver3(element_centre_x,element_centre_y,element_centre_z,fibre_U,fibre_V,fibre_W,10,'ShowArrowHead','off');
    %quiver3(element_centre_x,element_centre_y,element_centre_z,-fibre_U,-fibre_V,-fibre_W,10,'ShowArrowHead','off');
    quiver3(element_centre_x,element_centre_y,element_centre_z,-fibre_U,-fibre_V,-fibre_W,5,'ShowArrowHead','off',Color= fibre_colour,LineWidth=2);

     if PND == 1
         h=text(input.ND(node1,2) , input.ND(node1,3), input.ND(node1,4), num2str(node1)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node2,2) , input.ND(node2,3), input.ND(node2,4), num2str(node2)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node3,2) , input.ND(node3,3), input.ND(node3,4), num2str(node3)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node4,2) , input.ND(node4,3), input.ND(node4,4), num2str(node4)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node5,2) , input.ND(node5,3), input.ND(node5,4), num2str(node5)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node6,2) , input.ND(node6,3), input.ND(node6,4), num2str(node6)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node7,2) , input.ND(node7,3), input.ND(node7,4), num2str(node7)); set(h,'FontSize',8); 
         set(h,'Color','k');
         h=text(input.ND(node8,2) , input.ND(node8,3), input.ND(node8,4), num2str(node8)); set(h,'FontSize',8); 
         set(h,'Color','k');
             x1= input.ND(node1,2)/8  +  input.ND(node2,2)/8 +  input.ND(node3,2)/8  + input.ND(node4,2)/8 +...
                 input.ND(node5,2)/8  +  input.ND(node6,2)/8 +  input.ND(node7,2)/8  + input.ND(node8,2)/8;
             y1= input.ND(node1,3)/8  +  input.ND(node2,3)/8 +  input.ND(node3,3)/8  + input.ND(node4,3)/8 +...
                 input.ND(node5,3)/8  +  input.ND(node6,3)/8 +  input.ND(node7,3)/8  + input.ND(node8,3)/8;
             z1= input.ND(node1,4)/8  +  input.ND(node2,4)/8 +  input.ND(node3,4)/8  + input.ND(node4,4)/8 +...
                 input.ND(node5,4)/8  +  input.ND(node6,4)/8 +  input.ND(node7,4)/8  + input.ND(node8,4)/8;
             h=text(x1 , y1, z1, num2str(input.EL(i,1))); set(h,'FontSize',6); set(h,'Color','m');
     end;
 end;
 hold on;

 maxX = max(input.ND(:,2));  minX = min(input.ND(:,2));
 maxY = max(input.ND(:,3));  minY = min(input.ND(:,3));
 maxZ = max(input.ND(:,4));  minZ = min(input.ND(:,4));
 labx = (maxX / 9); laby = (maxY / 9); labz = (maxZ / 9);
 labx = min([labx laby labz]); laby = labx;  labz=labx;
%     labx = 0.3
%commented out lines that plot constraints and loads
% for i=1:size(input.LOAD)
%    node_i = find(input.ND(:,1)==input.LOAD(i,1));
%    if input.LOAD(i,2)~=0 & input.LOAD(i,2)>0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)+labx],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'r-','LineWidth',3); hold on;
%    end;
%    if input.LOAD(i,2)~=0 & input.LOAD(i,2)<0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)-labx],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'r-','LineWidth',3); hold on;
%    end;
%    if input.LOAD(i,3)~=0 & input.LOAD(i,3)>0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)+laby],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'r-','LineWidth',3); hold on;
%    end;
%    if input.LOAD(i,3)~=0 & input.LOAD(i,3)<0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)-laby],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'r-','LineWidth',3); hold on;
%    end;
%    if input.LOAD(i,4)~=0 & input.LOAD(i,4)>0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)+labz],'r-','LineWidth',3); hold on;
%    end;
%    if input.LOAD(i,4)~=0 & input.LOAD(i,4)<0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)-labz],'r-','LineWidth',3); hold on;
%    end;
% end;
%  labx = (maxX / 9); laby = (maxY / 9); labz = (maxZ / 9);
%  labx = min([labx laby labz]); laby = labx;  labz=labx;
% 
% for i=1:size(input.CON)
%    node_i = find(input.ND(:,1)==input.CON(i,1));
%    if input.CON(i,2)==0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)-labx],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'g-','LineWidth',2); hold on;
%    end;
%    if input.CON(i,3)==0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)-laby],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'g-','LineWidth',2); hold on;
%    end;
%    if input.CON(i,4)==0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)-labz],'g-','LineWidth',2); hold on;
%    end;
%    if input.CON(i,4)==0 %| input.CON(i,5)==0 %| input.CON(i,6)==0
%       plot3([input.ND(node_i,2) input.ND(node_i,2)],[input.ND(node_i,3) input.ND(node_i,3)],...
%          [input.ND(node_i,4) input.ND(node_i,4)],'gs','LineWidth',2); hold on;
%    end;
% end;

rotate3d(gca); set(gcf,'Pointer','arrow');
axis([(minX-labx) (maxX+labx) (minY-laby) (maxY+laby) (minZ-labz) (maxZ+labz)]); % plot scale
hold off    
end
