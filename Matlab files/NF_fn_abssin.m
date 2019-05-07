function [s_0,badj] = NF_fn_abssin(nnodes,thetas,alpha)
%% Node filtering model with edges determined by f(n) = sin(theta(n))

% Inputs:
    % nnodes        Number of nodes
    % thetas        2x1 vector of start and end theta, where theta =
    %                linspace(thetas(1),thetas(2),nnodes)
    % alpha         scaling factor
    
% Outputs:
    % s_0           Node ordering (1:nnodes)
    % badj          Binary adjacency matrix in s_0 order.
    
%
%
%% Main function:
    
    
t = linspace(thetas(1),thetas(2),nnodes);
% Preallocate
badj = zeros(nnodes);

% Run model
for n = 1:nnodes
    
    for j = 1:(n-1)
        
        % We will add an edge with probability alpha*abs(sin(t(n)))
        
        r = rand;
        
        if r < alpha*abs(sin(t(n)))
            
            badj(n,j) = 1;
            badj(j,n) = 1;
            
        end
        %disp(j)
    end
    
    
end
        
    
s_0 = 1:nnodes;

end

