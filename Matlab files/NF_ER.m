function [s_0,badj] = NF_constantProb(nnodes,p)
%% Node filtering model with edges determined by f(n) = p, a constant

% Inputs:
    % nnodes        Number of nodes
    
 
% Outputs:
    % s_0           Node ordering (1:nnodes)
    % badj          Binary adjacency matrix
    
    
% Preallocate
badj = zeros(nnodes);


% Run model
for n = 1:nnodes
    
    for j = 1:(n-1)
        
        % We will add an edge with probability 1/n^d
        
        r = rand;
        
        if r < p
            
            badj(n,j) = 1;
            badj(j,n) = 1;
            
        end
        %disp(j)
    end
    
    
end
        
    
s_0 = 1:nnodes;

end

