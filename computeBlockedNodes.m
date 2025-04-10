function B = computeBlockedNodes(G, centerNode, visualize)
% computeBlockedNodesRelaxed computes a static blocked node set for each node,
% using a relaxed criterion for allowed forwarding options, and optionally
% visualizes the resulting DAG subgraph.
%
%   B = computeBlockedNodesRelaxed(G, centerNode, visualize)
%
%   Inputs:
%       G         - A MATLAB graph object.
%       centerNode- The index of the center node.
%       visualize - (Optional) Boolean flag; if true, the function visualizes
%                   the DAG subgraph. Default is false.
%
%   Output:
%       B         - A cell array (length = number of nodes). For each node i,
%                   B{i} is a vector containing the indices of neighbors that
%                   are blocked for forwarding from node i.
%
%   Method:
%     1. Compute the shortest-path (hop-count) distances from each node to the center.
%     2. For each node i, define the allowed set A(i) as all neighbors j satisfying:
%            d(j) < d(i), OR
%            d(j) = d(i) and j < i  (using node indices as a tie-breaker)
%     3. Then, B(i) = Neighbors(i) \ A(i).
%     4. Optionally, if visualize is true, plot a directed graph where the edges
%        represent the allowed forwarding directions (i.e. the edges from i to j
%        for j in A(i)).
%
if nargin < 3
    visualize = false;
end

N = numnodes(G);

% Compute shortest-path distances (in hops) from each node to the center.
d = distances(G, centerNode);  % d(i) is the distance from node i to centerNode

B = cell(N, 1);
allowedEdges = [];  % Will store allowed edges as rows [i, j]

for i = 1:N
    nbrs = neighbors(G, i);
    allowed = [];
    for j = nbrs'
        if d(j) < d(i)
            allowed = [allowed; j];  %#ok<AGROW>
        elseif d(j) == d(i) && j < i
            allowed = [allowed; j];  %#ok<AGROW>
        end
    end
    % Blocked set: neighbors that are not allowed.
    B{i} = setdiff(nbrs, allowed);
    
    % For visualization, record allowed edges.
    for j = allowed'
        allowedEdges = [allowedEdges; i, j];  %#ok<AGROW>
    end
end

if visualize
    % Plot the original undirected graph.
    figure;
    h = plot(G, 'EdgeColor', [0.6,0.6,0.6], 'LineWidth', 1);
    hold on;
    title('Original Graph with Allowed Forwarding DAG Overlay');
    
    % Extract node coordinates (assumes G.Nodes has fields 'X' and 'Y')
    if isfield(G.Nodes, 'X') && isfield(G.Nodes, 'Y')
        X = G.Nodes.X;
        Y = G.Nodes.Y;
    else
        % If coordinates are not available, use a default layout.
        p = plot(G);
        X = p.XData;
        Y = p.YData;
        delete(p);
    end
    
    % Overlay allowed edges as bold red arrows.
    for idx = 1:size(allowedEdges,1)
        i = allowedEdges(idx,1);
        j = allowedEdges(idx,2);
        x1 = X(i); y1 = Y(i);
        x2 = X(j); y2 = Y(j);
        dx = x2 - x1;
        dy = y2 - y1;
        % Plot an arrow using quiver (scaling factor 0 so that dx,dy are taken as is).
        quiver(x1, y1, dx, dy, 0, 'LineWidth', 2, 'Color', 'r', 'MaxHeadSize', 0.5);
    end
    hold off;
end
end