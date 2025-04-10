function NetworkInfo = TopoGen(TopoType)
% TOPOGEN: Generates a graph according to the specified topology type.
%
%   [G, centerNode] = TopoGen(N_node_origin, TopoType, TopoPara)
%
%   Inputs:
%       N_node_origin - Approximate desired number of nodes.
%       TopoType      - Type of topology: 'grid' or 'er' (Erdos-Rényi).
%       TopoPara      - Parameter for the topology. For 'er', it is the
%                       link density (probability of link existence).
%
%   Outputs:
%       G         - Graph object representing the network.
%       centerNode- Node index with the minimum average shortest-path
%                   distance to all other nodes.
%
% Generate the graph based on the selected topology.
if strcmp(TopoType, 'grid')
    NetworkInfo = init_2d_grid();
elseif strcmp(TopoType, 'demo')
    [G, ~] = init_demo_topo();
elseif strcmp(TopoType, 'er')
    [G, ~] = init_ER_graph(N_node_origin, TopoPara);
else
    error('Unsupported topology type. Use ''grid'' or ''er''.');
end


G = NetworkInfo.G;
% Display the generated graph with node positions and highlight the center node.
figure;
p = plot(G, 'NodeLabel', {});
if isfield(G.Nodes, 'X') && isfield(G.Nodes, 'Y')
    p.XData = G.Nodes.X;
    p.YData = G.Nodes.Y;
end
title(sprintf('Topology: %s', TopoType));
hold on;
if isfield(G.Nodes, 'X') && isfield(G.Nodes, 'Y')
    plot(G.Nodes.X(centerNode), G.Nodes.Y(centerNode), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
hold off;
end

%% Function to generate a 2-D square grid graph with bi-directional links
function [G, centerNode] = init_demo_topo()
% a 3-node triangle, one server and 2 APs
coords = [0,sqrt(3); -1,0 ; 1,0];
edges = [1,2;2,1;1,3;3,1;2,3;3,2];
% Create the undirected graph. (graph() treats the input as undirected.)
G = graph(edges(:,1), edges(:,2));

% Add node coordinates as properties (useful for visualization).
G.Nodes.X = coords(:,1);
G.Nodes.Y = coords(:,2);

% Compute all-pairs shortest-path distances and select the node with
% minimum average distance as the center node.
D = distances(G);
avgDist = mean(D, 2);
[~, centerNode] = min(avgDist);

end

%% Function to generate a 2-D square grid graph with bi-directional links
function NetworkInfo = init_2d_grid()
% INIT_2D_GRID Generates a 2-D square grid graph with bi-directional links.
%
%   [G, centerNode] = init_2d_grid(N_node) returns an undirected graph G.
%   The grid will have total nodes equal to the square number closest to N_node.
%   Reverse edges are explicitly added so that each link is bi-directional.

% Determine grid size: choose the square number (m^2) closest to N_node.
N_node = 9;

% start
m0 = floor(sqrt(N_node));
m1 = ceil(sqrt(N_node));
if abs(m0^2 - N_node) <= abs(m1^2 - N_node)
    grid_size = m0;
else
    grid_size = m1;
end
N_node = grid_size^2;

% Generate other parameters
N_app = 5;                   % Number of applications
N_model = 2*ones(1,N_app);   % Number of models for each application (can vary)
maxN_model = max(N_model);   % Maximum number of models among apps
% Input rates for each application at each node
InputRate = 1 * ones(N_node, N_app);

LinkCostPara = zeros(N_node);  % Initialize full matrix
linkCapRef = 2*max(sum(InputRate,2)); % used to construct link capacities
CompCostPara = zeros(1, N_node);  % Computation capacities at nodes

ModelSize_ref = 1;
%ModelSize = ModelSize_ref * ones(N_app, maxN_model);  % Assume all models are of size 1
ModelSize = [ones(N_app,1),3*ones(N_app,1)];
eta = 1;                                % Tradeoff weight between delay and utility
%Reward = ones(N_app, maxN_model);       % Unit utility for all models
Reward = [0.1*ones(N_app,1),0.5*ones(N_app,1)];
L_req = 0.1 * ones(1, N_app);   % request size for each application
L_res  = 0.9 * ones(1, N_app);   % response size for each application
NodeCap = 3* ModelSize_ref * ones(1, N_node);   % storage size

% Create grid coordinates
[X, Y] = meshgrid(1:grid_size, 1:grid_size);
coords = [X(:), Y(:)];

% Build the list of directed edges: connect each node to its right and bottom neighbors.
edges = [];
for i = 1:grid_size
    for j = 1:grid_size
        current = (i-1)*grid_size + j;

        layer = min(i,j);
        % Connect to the right neighbor, if it exists.
        if j < grid_size
            right = current + 1;
            edges = [edges; current, right];  %#ok<AGROW>

            LinkCostPara(current,right) = linkCapRef * layer;
        end
        % Connect to the bottom neighbor, if it exists.
        if i < grid_size
            down = current + grid_size;
            edges = [edges; current, down];  %#ok<AGROW>

            LinkCostPara(current,down) = linkCapRef * layer;
        end
        CompCostPara(current) = linkCapRef * layer;
    end
end

% Now, add reverse edges so that each link is explicitly bi-directional.
reverseEdges = [edges(:,2), edges(:,1)];
allEdges = [edges; reverseEdges];
% Remove duplicate rows (if any)
allEdges = unique(allEdges, 'rows');

% Create the undirected graph. (graph() treats the input as undirected.)
G = graph(allEdges(:,1), allEdges(:,2));

% Add node coordinates as properties (useful for visualization).
G.Nodes.X = coords(:,1);
G.Nodes.Y = coords(:,2);

% Compute all-pairs shortest-path distances and select the node with
% minimum average distance as the center node.
centerNode = current;

CompCostPara(centerNode) = 1e6;
NodeCap(centerNode) = 1e6;  % Set center node capacity very high

LinkCostPara = LinkCostPara + LinkCostPara.';

NetworkInfo.G = G;
NetworkInfo.centerNode = centerNode;
NetworkInfo.N_app = N_app;
NetworkInfo.N_model = N_model;
NetworkInfo.InputRate = InputRate;
NetworkInfo.LinkCostPara = LinkCostPara;
NetworkInfo.CompCostPara = CompCostPara;
NetworkInfo.ModelSize = ModelSize;
NetworkInfo.Reward = Reward;
NetworkInfo.L_req = L_req;
NetworkInfo.L_res = L_res;
NetworkInfo.NodeCap = NodeCap;
NetworkInfo.eta = eta;

end

%% Function to generate an Erdos-Rényi graph ensuring connectivity
function [G, centerNode] = init_ER_graph(N_node, linkProb)
% INIT_ER_GRAPH Generates a connected Erdos-Rényi graph.
%
%   [G, centerNode] = init_ER_graph(N_node, linkProb) returns an undirected graph G
%   where each possible edge exists with probability linkProb.
%   Nodes are assigned random coordinates in the unit square [0,1]x[0,1].
%   The function ensures that the generated graph is connected.
%
% Note: The Erdos-Rényi graph is already made symmetric.

maxTries = 1000;
isConnected = false;
count = 0;

while ~isConnected && count < maxTries
    count = count + 1;
    % Generate an upper-triangular random matrix and symmetrize it.
    A = triu(rand(N_node) < linkProb, 1);
    A = A + A';

    % Create the graph.
    G = graph(A);

    % Check connectivity using connected components.
    comp = conncomp(G);
    if numel(unique(comp)) == 1
        isConnected = true;
    end
end

if ~isConnected
    error('Failed to generate a connected Erdos-Rényi graph after %d attempts.', maxTries);
end

% Assign random coordinates to nodes in the unit square.
coords = rand(N_node, 2);
G.Nodes.X = coords(:,1);
G.Nodes.Y = coords(:,2);

% Compute all-pairs shortest-path distances and select the node with
% minimum average distance as the center node.
D = distances(G);
avgDist = mean(D, 2);
[~, centerNode] = min(avgDist);
end
