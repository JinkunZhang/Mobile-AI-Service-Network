function [s_new, Phi_new] = updateRoutingAndSelectionLP(G, y, InputRate, LinkCostPara, Reward)
% updateRoutingAndSelectionLP assigns routing and model selection using
% shortest path (zero-flow delay) to nearest model-hosting node.

N_node = size(InputRate,1);
N_app = size(InputRate,2);
maxN_model = size(y,3);

% Initialize
s_new = zeros(size(y));
Phi_new = zeros(N_node, N_app, maxN_model, N_node);

% Build graph with static delay weights
delay_graph = digraph(LinkCostPara);

for k = 1:N_app
    for i = 1:N_node
        % Skip if no input at this node for this app
        if InputRate(i,k) == 0
            continue;
        end

        % For each model m, find the nearest node that hosts it
        best_dist = Inf;
        best_m = 1;
        best_path = [];
        
        for m = 1:maxN_model
            % Find nodes hosting (k,m)
            host_nodes = find(y(:,k,m) > 0);
            if isempty(host_nodes), continue; end
            
            % Compute shortest path from i to any host node
            min_dist = Inf;
            best_host = -1;
            best_local_path = [];

            for j = host_nodes(:)'
                [path, d] = shortestpath(delay_graph, i, j);
                if d < min_dist
                    min_dist = d;
                    best_host = j;
                    best_local_path = path;
                end
            end
            
            if min_dist < best_dist
                best_dist = min_dist;
                best_m = m;
                best_path = best_local_path;
            end
        end

        % Assign selection and routing for best model m
        s_new(i,k,best_m) = 1;
        if length(best_path) >= 2
            for idx = 1:length(best_path)-1
                j_from = best_path(idx);
                j_to = best_path(idx+1);
                Phi_new(i,k,best_m,j_to) = 1;  % deterministic
            end
        end
    end
end
end
