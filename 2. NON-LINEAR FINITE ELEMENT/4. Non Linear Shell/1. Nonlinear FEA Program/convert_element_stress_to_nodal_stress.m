function [storeF, storeStress] = convert_element_stress_to_nodal_stress(F, Stress, input, load_step)

% element_stress gives elements at each local node (and centre of element)
% for each element for a given load step (only one load step should be fed
% in here)

num_nodes = numel(input.ND(:,1));
num_eles  = numel(input.EL(:,1));

element_F = F(:,:,load_step);
element_S = Stress(:,:,load_step);

nGradientComponents = 9;
nodesPerElem = 8;

storeF = zeros(num_nodes,nGradientComponents);
storeStress = zeros(num_nodes,nGradientComponents);

connectedElements = zeros(nodesPerElem,2);

% Perform conversion to principal material coodinates...
convert_to_FSN = 1;

if convert_to_FSN == 1
    element_F_new = element_F*0;
    element_S_new = element_S*0;
    for elem = 1:num_eles
        f = input.FIBRES.f(elem,:);
        s = input.FIBRES.s(elem,:);
        n = input.FIBRES.n(elem,:);
        Q = [f; s; n]';
        for node = 1:nodesPerElem
            Fvec = element_F(elem,(node-1)*nGradientComponents+1:node*nGradientComponents);
            Fmat = V2M(Fvec');
            F_fsn = Q.'*Fmat*Q;

            S_fsn = Get_S_CC(F_fsn,input.model,input.FIBRES);
            element_F_new(elem,(node-1)*nGradientComponents+1:node*nGradientComponents) = M2V(F_fsn)';
            element_S_new(elem,(node-1)*nGradientComponents+1:node*nGradientComponents) = M2V(S_fsn)';
        end
    end
    
    element_F = element_F_new;
    element_S = element_S_new;

end





for node = 1 : num_nodes
    
    % This could be improved by using another loop and arrays..
    i = 0;
    for localNode = 1 : nodesPerElem
        iElement = find(input.EL(:,2+localNode) == node);
        if isempty(iElement) ~= 1
            i = i + 1;
            connectedElements(i,:) = [iElement, localNode];
        end
    end
    numConnected = i;

    % Get F components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = element_F(element,nGradientComponents*(localNode-1)+component);
        end
        storeF(node,component) = mean(values);
    end

    % Get tau components
    values = zeros(numConnected,1);
    for component = 1 : nGradientComponents
        for element = 1 : numConnected
            localNode = connectedElements(element,2);
            values(element) = element_S(element,nGradientComponents*(localNode-1)+component);
        end
        storeStress(node,component) = mean(values);
    end
           
end    