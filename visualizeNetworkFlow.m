function visualizeNetworkFlow(G, F_total,MaxFlow)
% visualizeNetworkFlow plots the network topology with edge line widths
% proportional to the flow in F_total.
%
%   visualizeNetworkFlow(G, F_total)
%
%   Inputs:
%       G       - MATLAB graph object.
%       F_total - [N x N] matrix of aggregated flows on each link.
%
% This function uses the node coordinates stored in G.Nodes.X and G.Nodes.Y.
% It loops through each edge in the graph and draws a red line with width
% proportional to F_total(i,j) for edge (i,j).

% Get node coordinates
if isfield(G.Nodes, 'X') && isfield(G.Nodes, 'Y')
    X = G.Nodes.X;
    Y = G.Nodes.Y;
else
    % If coordinates are not available, use a default layout.
    p_temp = plot(G);
    X = p_temp.XData;
    Y = p_temp.YData;
    close(gcf);
end

% Define a scaling factor for line width (adjust as needed)
scalingFactor = 20 / MaxFlow;

% Create a new figure (or hold on to current axes)
figure;
hold on;

% Plot nodes as black circles.
plot(X, Y, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Loop through each edge in the graph.
EdgeTable = G.Edges;
for e = 1:height(EdgeTable)
    i = EdgeTable.EndNodes(e,1);
    j = EdgeTable.EndNodes(e,2);
    flow = F_total(i,j);
    % Compute line width; clip to a minimum of 1 and a maximum of 20.
    lw = scalingFactor * flow;
    lw = max(1, min(lw, 20));
    % Draw the edge as a red line.
    plot([X(i), X(j)], [Y(i), Y(j)], 'r-', 'LineWidth', lw);
end

title('Network Flow Visualization');
axis equal;
hold off;
end
