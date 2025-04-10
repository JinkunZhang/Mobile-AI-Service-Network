function y = y_init_gen(N_node, N_app, N_model, maxN_model, NodeCap, init_type)
% Y_INIT_GEN generates a valid initial model placement (y) for each node.
%
%   y = y_init_gen(N_node, N_app, N_model, maxN_model, NodeCap, init_type)
%
%   Inputs:
%       N_node   - Number of nodes in the network.
%       N_app    - Number of applications.
%       N_model  - A vector of length N_app specifying the number of models for each application.
%       maxN_model - Maximum number of models across all applications.
%       NodeCap  - Vector of node capacities (1 x N_node). For non-center nodes,
%                  the full capacity will be allocated to models.
%       init_type- Initialization type: currently, use 'uniform' for uniform distribution.
%
%   Output:
%       y        - A 3D array of size [N_node x N_app x maxN_model].
%
%   For 'uniform' initialization:
%       - For each non-center node i, the total capacity NodeCap(i) is split equally
%         among all models (across all applications). That is, for each valid (k,m):
%              y(i,k,m) = NodeCap(i) / (sum(N_model))
%       - For the center node (global CENTER_NODE), y is set to 1 for each valid model,
%         meaning that each model is permanently stored.
%   For 'cloud':
%       - set all y=0 except central node where all y = 1

    % Set default initialization type if not provided.
    if nargin < 6
        init_type = 'uniform';
    end

    % Global variable for center node index must be set externally.
    global CENTER_NODE;
    if isempty(CENTER_NODE)
        error('Global variable CENTER_NODE is not defined.');
    end

    y = zeros(N_node, N_app, maxN_model);
    
    switch lower(init_type)
        case 'uniform'
            % Compute total number of valid models across all applications.
            totalModels = sum(N_model);
            for i = 1:N_node
                for k = 1:N_app
                    for m = 1:maxN_model
                        if m <= N_model(k)
                            if i == CENTER_NODE
                                % For the center node, store every model permanently.
                                y(i, k, m) = 1;
                            else
                                % For non-center nodes, use full capacity,
                                % equally split among all valid models.
                                y(i, k, m) = NodeCap(i) / totalModels;
                            end
                        else
                            y(i, k, m) = 0;
                        end
                    end
                end
            end
        case 'cloud'
            % Compute total number of valid models across all applications.
            totalModels = sum(N_model);
            for i = 1:N_node
                for k = 1:N_app
                    for m = 1:maxN_model
                        if m <= N_model(k)
                            if i == CENTER_NODE
                                % For the center node, store every model permanently.
                                y(i, k, m) = 1;
                            else
                                y(i, k, m) = 0;
                            end
                        else
                            y(i, k, m) = 0;
                        end
                    end
                end
            end
        otherwise
            error('Unknown initialization type: %s', init_type);
    end

    % Validation: For each non-center node, ensure the total allocated equals NodeCap(i).
    for i = 1:N_node
        if i ~= CENTER_NODE
            total_y = sum(sum(y(i, :, :)));
            if (total_y - NodeCap(i)) > 1e-6
                error('Initial y exceeds capacity constraint at node %d: total = %.4f, capacity = %.2f', i, total_y, NodeCap(i));
            end
        end
    end

end