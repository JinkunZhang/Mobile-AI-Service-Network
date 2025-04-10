function [s_init, phi_init] = s_phi_init_gen(G, centerNode, N_app, N_model, maxN_model, y_init, initOption)
% S_PHI_INIT_GEN generates initial values for s and phi.
%
% s_init: [N_node x N_app x maxN_model], where for each node i and app k,
%         s_init(i,k,1:N_model(k)) = 1/N_model(k) and the rest are 0.
%
% phi_init: [N_node x N_app x maxN_model x N_node]. For each node i (≠ center):
%   - If initOption is 'shortestPath' (default):
%         Let nextHop be the second node on the shortest path from i to center.
%         Then, for each app k and model m, set:
%            phi_init(i,k,m,nextHop) = 1 - y_init(i,k,m)
%   - If initOption is 'even':
%         Let allowed = neighbors(i) \ B{i}, where B = computeBlockedNodes(G, centerNode, false).
%         Then, for each app k and model m, set:
%            phi_init(i,k,m,allowed) = (1 - y_init(i,k,m)) / numel(allowed)
%   - If initOption is 'random':
%         Let allowed = neighbors(i) \ B{i} (or all neighbors if empty).
%         Then, for each app k and model m, generate a random vector (with a fixed seed)
%         of length numel(allowed) and set:
%            phi_init(i,k,m,allowed) = (1 - y_init(i,k,m)) * (random vector normalized to sum to 1)
%
% For the center node, phi is all zeros.
%
% If initOption is not provided, it defaults to 'shortestPath'.

if nargin < 7 || isempty(initOption)
    initOption = 'shortestPath';
end

N_node = numnodes(G);
s_init = zeros(N_node, N_app, maxN_model);
phi_init = zeros(N_node, N_app, maxN_model, N_node);

% Initialize s_init: equally distribute among valid models.
for i = 1:N_node
    for k = 1:N_app
        for m = 1:maxN_model
            if m <= N_model(k)
                s_init(i,k,m) = 1 / N_model(k);
            else
                s_init(i,k,m) = 0;
            end
        end
    end
end

switch lower(initOption)
    case 'shortestpath'
        % For each node (except center), use the shortest path to center.
        for i = 1:N_node
            if i == centerNode
                phi_init(i,:,:,:) = 0;
            else
                path = shortestpath(G, i, centerNode);
                if length(path) >= 2
                    nextHop = path(2);
                else
                    nextHop = i;
                end
                for k = 1:N_app
                    for m = 1:maxN_model
                        phi_init(i,k,m,nextHop) = 1 - y_init(i,k,m);
                    end
                end
            end
        end
        
    case 'even'
        % Compute blocked node sets for each node.
        BlockNodes = computeBlockedNodes(G, centerNode, false);
        % For each node (except center), distribute load evenly among allowed neighbors.
        for i = 1:N_node
            if i == centerNode
                phi_init(i,:,:,:) = 0;
            else
                nbrs = neighbors(G, i);
                allowed = setdiff(nbrs, BlockNodes{i});
                if isempty(allowed)
                    allowed = nbrs;
                end
                for k = 1:N_app
                    for m = 1:maxN_model
                        phi_init(i,k,m,allowed) = (1 - y_init(i,k,m)) / numel(allowed);
                    end
                end
            end
        end
        
    case 'random'
        % Compute blocked node sets for each node.
        BlockNodes = computeBlockedNodes(G, centerNode, false);
        % Set a fixed random seed for reproducibility.
        rng(123456);
        for i = 1:N_node
            if i == centerNode
                phi_init(i,:,:,:) = 0;
            else
                nbrs = neighbors(G, i);
                allowed = setdiff(nbrs, BlockNodes{i});
                if isempty(allowed)
                    allowed = nbrs;
                end

                 % Generate a random vector for allowed neighbors.
                randVec = rand(length(allowed), 1);
                randVec = randVec / sum(randVec); % Normalize to sum to 1.

                for k = 1:N_app
                    for m = 1:maxN_model
                       
                        % Assign to phi_init.
                        phi_init(i,k,m,allowed) = (1 - y_init(i,k,m)) * randVec;
                    end
                end
            end
        end
        
    otherwise
        error('Unknown initialization option: %s. Use ''shortestPath'', ''even'', or ''random''.', initOption);
end
end
