function [s_0,badj] = NF_proportionalProb(nNodes,d)
%% Node filtering model with edges determined by f(n) = (n/nNodes)^d

% Inputs:
    % nNodes        Number of nodes.
    % d             Expondent, require d>0.
    
% Outputs:
    % s_0           Node ordering (1:nnodes).
    % badj          Binary adjacency matrix in s_0 order.
    
%
% Reference: Ann E. Sizemore, Elisabeth A. Karuza, Chad Giusti, Danielle S.
% Bassett, Knowledge gaps in the early growth of semantic networks.
% Submitted.
%
%% Main function:
    
% Preallocate
badj = zeros(nNodes);

if d<0
    error('d should be at least 0')
end

% Run model
for n = 2:nNodes
    
    for j = 1:(n-1)
        
        % We will add an edge with probability (n/nNodes)^d
        
        r = rand;
        
        if r < (n/nNodes)^d
            
            badj(n,j) = 1;
            badj(j,n) = 1;
        end
        
    end
    
    
end
        
    
s_0 = 1:nNodes;

end

