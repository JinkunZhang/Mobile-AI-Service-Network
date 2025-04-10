function delta_cell = computeDelta(Phi, dJdF, y, G_workload, CompCostType, CompCostPara, LinkPara,L_req,L_res)
% computeDelta computes the exact delay sensitivity δ for each node, application, 
% and model using the outgoing routing probabilities and the derivative of J 
% with respect to the outgoing flow, dJdF.
%
%   delta_cell = computeDelta(Phi, dJdF, y, G_workload, CompCostType, CompCostPara, LinkPara)
%
%   Inputs:
%       Phi         - [N_node x N_app x maxN_model x N_node] routing probabilities.
%                     Phi(i,k,m,j) is the fraction of flow forwarded from node i to node j.
%       dJdF      - [N_node x N_node] matrix representing ∂J/∂F^o.
%       y           - [N_node x N_app x maxN_model] local processing fractions.
%       G_workload  - [N_node x 1] vector of node workloads (used for computing the 
%                     computation derivative).
%       CompCostType- String specifying the computation cost type ('linear','queue','taylor').
%       CompCostPara- Computation capacity parameter (vector of length N_node).
%       LinkPara    - [N_node x N_node] matrix of link parameters (nonzero indicates a valid link).
%
%   Output:
%       delta_cell  - Cell array of size [N_app x maxN_model], where each cell is an
%                     [N_node x 1] vector representing δ for that (k,m), computed as:
%
%         (I - Phi_out) * δ = δ_local + b,
%
%         where δ_local(i) = y(i,k,m) * CompDelayDerivative(G_workload(i),...)
%         and b(i) = sum_{j in O(i)} Phi(i,k,m,j) * dJdF(i,j).
%
%   This exactly mimics the messaging protocol described in the paper.
%

N_node = size(Phi, 1);
N_app = size(Phi, 2);
maxN_model = size(Phi, 3);

% Compute computation derivative at each node.
comp_deriv = CompDelayDerivative(G_workload, CompCostType, CompCostPara);

delta_cell = cell(N_app, maxN_model);
I = eye(N_node);

for k = 1:N_app
    for m = 1:maxN_model
        % Extract the routing matrix for branch (k,m):
        % Phi_km(i,j) = Phi(i,k,m,j)
        Phi_km = squeeze(Phi(:, k, m, :));  % [N_node x N_node]

        Req_size = L_req(k);
        Res_size = L_res(k);
        
        % Compute local delay sensitivity δ_local for each node:
        % δ_local(i) = y(i,k,m) * comp_deriv(i)
        delta_local = squeeze(y(:, k, m)) .* comp_deriv;
        
        % Compute b vector for each node i:
        % b(i) = sum_{j in O(i)} Phi(i,k,m,j) * dJdF(i,j)
        % We determine valid neighbors as those with LinkPara(i,j) ~= 0.
        b = zeros(N_node, 1);
        for i = 1:N_node
            validNbrs = find(LinkPara(i,:) ~= 0);
            b(i) = sum(Phi_km(i, validNbrs) .* (Req_size * dJdF(i, validNbrs) + Res_size * dJdF(validNbrs,i).'));
        end
        
        % Solve the linear system: (I - Phi_km) * δ = δ_local + b.
        rhs = delta_local + b;
        delta_cell{k, m} = (I - Phi_km) \ rhs;
    end
end
end
