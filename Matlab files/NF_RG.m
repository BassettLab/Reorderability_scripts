function [s_0,badj] = NF_RG(nnodes,ep,dim)
%% Node filtering model with edges determined by f(n) = 1/n^d

% Inputs:
    % nnodes        Number of nodes
    % ep            network threshold
    % dim           dimension of unit cube
    
% Outputs:
    % s_0           Node ordering (1:nnodes)
    % badj          Binary adjacency matrix
    


points = rand(nnodes,dim);

% filter through x axis
[~,inds1] = sort(points(:,1),'ascend');
points_sorted = points(inds1,:);

% Now ordering can be 1:nNodes
D = pdist(points_sorted);
adj = 1./squareform(D);
adj(logical(eye(nnodes))) = 0;

badj = thresholdMatDensity(adj,ep);
badj(badj>0) = 1;

s_0 = 1:nnodes;
    

end

