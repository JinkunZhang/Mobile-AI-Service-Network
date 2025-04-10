function G_workload = ComputeNodeWorkload(InputRate, s, f_req_cell, ModelSize)
% ComputeNodeWorkloadUpdated calculates the computation workload at each node
% using the exogenous request rates and the forwarded (incoming) flows.
%
%   G_workload = ComputeNodeWorkloadUpdated(InputRate, s, f_req_cell, ModelSize)
%
%   Inputs:
%       InputRate  - [N_node x N_app] matrix of input rates at each node.
%       s          - [N_node x N_app x maxN_model] array with model selection fractions.
%       f_req_cell - Cell array of size [N_app x maxN_model] from CalcRequestFlow.
%                    Each cell is an [N_node x N_node] matrix, where the (l,i) entry is 
%                    the flow from node l to node i for a given (k,m).
%       ModelSize  - [N_app x maxN_model] matrix where ModelSize(k,m) is the 
%                    computational cost per request \(\omega_{k,m}\) for model m of app k.
%
%   Output:
%       G_workload - [N_node x 1] vector of computation workloads at each node.
%
%   For each node i and each (k,m):
%       t_i^{k,m} = InputRate(i,k)*s(i,k,m) + sum_{l=1}^{N_node} f_req_cell{k,m}(l,i)
%   and then:
%       G(i) = \sum_{k=1}^{N_app}\sum_{m=1}^{maxN_model} ModelSize(k,m) * t_i^{k,m}.
%

N_node = size(InputRate, 1);
N_app = size(InputRate, 2);
maxN_model = size(s, 3);

G_workload = zeros(N_node, 1);

for k = 1:N_app
    for m = 1:maxN_model
        % Compute exogenous traffic r_{i,k,m}
        r_km = InputRate(:, k) .* squeeze(s(:, k, m));  % [N_node x 1]
        % Compute the forwarded (incoming) traffic at each node:
        % f_req_cell{k,m} is an [N_node x N_node] matrix. Its (l,i) entry is the flow from node l to node i.
        forwarded = sum(f_req_cell{k,m}, 1)';  % Sum over rows gives a [N_node x 1] vector.
        % Total effective traffic at each node for (k,m)
        t_km = r_km + forwarded;
        % Accumulate workload weighted by model cost.
        G_workload = G_workload + ModelSize(k, m) * t_km;
    end
end

end
