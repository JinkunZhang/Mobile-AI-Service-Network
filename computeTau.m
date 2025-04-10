function tau = computeTau(D0_cell, Lambda_Node, q_mat, D_prime)
% computeTau computes the tunneling-related delay derivative tau for each node,
% application, and model.
%
%   tau = computeTau(D0_cell, Lambda_Node, q_mat, d_prime)
%
%   Inputs:
%       D0_cell    - Cell array of size [N_app x maxN_model], where each cell is an
%                    [N_node x 1] vector of expected end-to-end delays D0 for that (k,m).
%       Lambda_Node- [1 x N_node] vector of node transition rates.
%       q_mat      - [N_node x N_node] matrix of movement probabilities (q_{ij}).
%       d_prime    - [N_node x N_node] matrix of the derivatives of one-way link delays (d'_{ij}).
%
%   Output:
%       tau        - [N_node x N_app x maxN_model] array, where for each node i, app k, and model m:
%
%              tau(i,k,m) = sum_{j=1}^{N_node} d_prime(i,j)*q_mat(i,j) * (1 - exp(-Lambda_Node(i)*D0_cell{k,m}(i))).
%
%   This term represents the marginal latency due to tunneling of exogenously
%   arrived requests at node i.
%
N_node = size(q_mat, 1);
[N_app, maxN_model] = size(D0_cell);

tau = zeros(N_node, N_app, maxN_model);

for k = 1:N_app
    for m = 1:maxN_model
        % Retrieve the expected delay vector for application k, model m.
        D0_vec = D0_cell{k, m};  % [N_node x 1]
        for i = 1:N_node
            % Compute the factor 1 - exp(-Lambda_i * D0(i,k,m))
            factor = 1 - exp(-Lambda_Node(i) * D0_vec(i));
            % Sum over all j (neighbors): d_prime(i,j) * q_mat(i,j)
            tau(i,k,m) = sum(D_prime(i,:) .* q_mat(i,:), 'omitnan') * factor;
        end
    end
end

end
