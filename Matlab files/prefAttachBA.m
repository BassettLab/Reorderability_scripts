function [ badj, s_0 ] = prefAttachBA( nNodes, m_0,m)
% Classic preferential attachment model (Barabasi-Albert) 
% Revised 09.18
%
%
% Inputs:
%       nNodes = number of nodes in output
%       M = number of connections to add at each timestep
%
% Output:
%       adj = binary matrix of output network
%
%
%
%%   Main function:

badj = zeros(nNodes);

% begin - random network
isConnected = 0;

while isConnected ~= 1
    randNet = rand(m_0)-eye(m_0);
    randNet = (randNet+randNet')/2;
    randNet(randNet<0.5)= 0;
    randNet(randNet>0) = 1;

    % ensure network is connected
    comps = get_components(randNet);
    isConnected = max(comps);
end

badj(1:m_0,1:m_0) = randNet;


%% add nodes 

for n = (m_0+1):nNodes
    
    P = sum(badj)./sum(badj(:));
    cP = cumsum(P);
    
    for m_i = 1:m

        % choose node i
        r = rand;
        node_i = find(r<cP,1,'first');
        
        while ismember(node_i,find(badj(n,:)))
            r = rand;
            cP = cumsum(P);
            node_i = find(r<cP,1,'first');
        end

        badj(n,node_i) = 1;
        badj(node_i,n) = 1;
    end


end

    
s_0 = 1:nNodes;


end

