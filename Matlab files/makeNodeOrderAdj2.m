function [jadj,adj,adjn] = makeNodeOrderAdj2(s,badj)
%% Make weighted matrix from node order and binary graph

% Given weighted nodes and a binary graph of edges between these nodes, we 
% will construct a weighted graph with the ordering of the edges as they
% would appear if we filtered by decreasing node weight. 

% INPUTS
    % s: a 1 x nnodes vector indicating node order
    % badj: binary adjacency matrix 
            
            
% OUTPUT
    % adj: weighted adjacency matrix with larger weights corresponding to
            % earlier-appearing edges (ordering induced by node weight ordering)
    % adjn: adj + noise (can change at lines 50-55)
    % svec_out: 1xnEdges vector noting which edge belongs to which node.

badj_reordered = badj(s,s);

nNodes = length(s);

valMat = ones(nNodes);

for col = 1:nNodes
    valMat(1:col,col) = valMat(1:col,col)*col;
    valMat(col,1:col) = valMat(col,1:col)*col;
end

jadj = badj_reordered.*valMat;
adj = 1./jadj;

jadj(jadj==0) = 2*nNodes;
jadj(logical(eye(nNodes))) = 0;

adj(isinf(adj)) = 0;

adjn = adj;
for i = 1:nNodes
    for j = i+1:nNodes
        adjn(i,j) = adjn(i,j) + 0.00001*rand;
        adjn(j,i) = adjn(i,j);
    end
end

end

