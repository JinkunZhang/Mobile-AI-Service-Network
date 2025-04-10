function hopDist = computeHopDistanceToCenter(G, centerNode)
% computeHopDistanceToCenter computes the hop distance from each node to a specified center node.
%
%   hopDist = computeHopDistanceToCenter(G, centerNode)
%
%   Inputs:
%       G           - MATLAB graph object
%       centerNode  - Scalar index of the central node (e.g., the server node)
%
%   Output:
%       hopDist     - [N_node x 1] vector, where hopDist(i) is the minimum
%                     number of hops from node i to centerNode

N_node = numnodes(G);
hopDist = inf(N_node, 1);

% Use MATLAB's built-in shortestpath graph function
for i = 1:N_node
    if i == centerNode
        hopDist(i) = 0;
    else
        [~, hopDist(i)] = shortestpath(G, i, centerNode);
    end
end
end
