%%

function [K0, DeltaP, initialLoad, finalLoad, maxInc] = pre_processing_NFEA(input, elementType, Lx, Ly, ~, CCASE, FCASE, totalDOFs, nodesPerElem, dofPerNode, numLoadSteps, finalStress)




%% Display details of defined elasticity problem
fprintf('\n   Problem Generated.');
fprintf('\n');
fprintf('\n   Summary of Problem:');
fprintf('\n      Body width  = %.2f mm',Lx);
fprintf('\n      Body height = %.2f mm',Ly);
fprintf('\n      Constriants:');
if CCASE == 1
    fprintf('\n         Bottom edge pinned.');
elseif CCASE == 2
    fprintf('\n         Centre node pinned.');
end
fprintf('\n      Loading:');
if FCASE == 1
    fprintf('\n         Top edge uniform load = %.2f N/mm2',finalStress);
elseif FCASE == 2
    fprintf('\n         Top edge uniform load    =  %.2f N/mm2',finalStress);
    fprintf('\n         Bottom edge uniform load = -%.2f N/mm2',finalStress);
end
fprintf('\n');
fprintf('\n   Summary of Discretisation:');
if elementType == 1
    fprintf('\n      2D Bilinear Quadrilateral Elements');
elseif elementType == 2
    fprintf('\n      3D Hexahedral Elements');
end
fprintf('\n      Number of Nodes:          %i',size(input.ND,1));
fprintf('\n      Number of Elements:       %i',size(input.EL,1));
fprintf('\n      Total Degrees of Freedom: %i',totalDOFs);
fprintf('\n');
fprintf('\n');


%%                           Plot of Discretisation

% This section produces a plot of the discretised mesh and applied boundary
% condtitions (appropriated from open-source code)

fig = figure;
set(fig,'Name','FEM Problem Mesh','NumberTitle','off');
plotDiscretisation(input,1,elementType)

%%                         Determine Initial Gradient

% This section determines the initial gradient of the elasticity problem.
% This is equivalent to the tangent stiffness matrix of each element 
% evaluated in the region of infinitesimal displacement gradients

fprintf('\n   Determining Initial Stiffness Matrix...');

K0 = Assemble_Initial_Stiffness(input,elementType,nodesPerElem,dofPerNode);

fprintf('\n   Completed.');
fprintf('\n');
fprintf('\n   Summary of Initial Stiffness Matrix:');
fprintf('\n      Total size:             %i x %i',dofPerNode*size(input.ND,1),dofPerNode*size(input.ND,1));
fprintf('\n      Total entries:          %i',(dofPerNode*size(input.ND,1))^2);
fprintf('\n      Total non-zero entries: %i',nnz(K0));
fprintf('\n');


%%                   Define Incremental Load Vector

% This section generates the final load vector based on input conditions,
% before using this to determine the global incremental load vector

fprintf('\n   Constructing incremental load vector...');

% Number of load steps
fprintf('\n      Total load steps:     %i',numLoadSteps);

% Final Load Vector
finalLoad = Final_Load_Vector(input);
fprintf('\n      Final load vector constructed');
fprintf('\n      Maximum Load:         %.3f N',max(max(input.LOAD(:,2:end))));

% Initial Load Vector
initialLoad = 0.*finalLoad;

% Incremental load vector
DeltaP = (1/numLoadSteps).*(finalLoad - initialLoad);
fprintf('\n      Incremental load vector determined');
fprintf('\n');
maxInc = max(max(input.LOAD(:,2:end)))/numLoadSteps;


end