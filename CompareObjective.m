function [J_total, AvgDelay, AvgReward] = CompareObjective(F_req, F_res, F_tun, InputRate, s, y, Reward, eta, CompCostType, CompCostPara, LinkCostType, LinkCostPara)
% CompareObjective compares the overall objective function computed with total flows 
% (F_total = F_req + F_res + F_tun) versus flows without tunneling (F_no_tun = F_req + F_res).
%
%   [J_total, J_no_tun, diffJ] = CompareObjective(F_req, F_res, F_tun, InputRate, s, y, Reward, eta, ...
%       CompCostType, CompCostPara, LinkCostType, LinkCostPara)
%
%   Inputs:
%       F_req         - [N_node x N_node] aggregated request flow matrix.
%       F_res         - [N_node x N_node] aggregated response flow matrix.
%       F_tun         - [N_node x N_node] aggregated tunneling flow matrix.
%       InputRate     - [N_node x N_app] matrix of input rates.
%       s             - [N_node x N_app x maxN_model] array of model selection fractions.
%       y             - [N_node x N_app x maxN_model] array of local processing fractions.
%       Reward        - [N_app x maxN_model] matrix of unit utilities.
%       eta           - Tradeoff weight between delay and utility.
%       CompCostType  - String specifying computation cost type ('linear','queue','taylor').
%       CompCostPara  - Vector of computation capacities (length = N_node).
%       LinkCostType  - String specifying link cost type ('linear','queue','taylor').
%       LinkCostPara  - [N_node x N_node] matrix of link capacities.
%
%   Outputs:
%       J_total  - Overall objective computed with tunneling (total flows).
%       J_no_tun - Overall objective computed without tunneling.
%       diffJ   - Difference: J_no_tun - J_total (how much extra cost/delay is introduced by tunneling).
%
%   The objective is defined as:
%       J = [eta*(total utility) - (total link cost) - (total computation cost)] / (total request rate)
%

N_node = size(F_req,1);
N_app = size(InputRate,2);
maxN_model = size(s,3);

%% Compute Total Flows and Link Delays
F_total = F_req + F_res + F_tun;
F_no_tun = F_req + F_res;

d_total = LinkDelay(F_total, LinkCostType, LinkCostPara);
d_no_tun = LinkDelay(F_no_tun, LinkCostType, LinkCostPara);

% Create a valid mask: only include links that are valid.
valid_mask = (LinkCostPara > 0) & ~isnan(d_total); %& ~isinf(d_total);
link_cost_total = sum(sum(F_total(valid_mask) .* d_total(valid_mask)));

valid_mask_no_tun = (LinkCostPara > 0) & ~isnan(d_no_tun);% & ~isinf(d_no_tun);
link_cost_no_tun = sum(sum(F_no_tun(valid_mask_no_tun) .* d_no_tun(valid_mask_no_tun)));

%% Compute Node Workload and Computation Cost
G_workload = zeros(N_node,1);
for i = 1:N_node
    for k = 1:N_app
        for m = 1:maxN_model
            G_workload(i) = G_workload(i) + InputRate(i,k) * s(i,k,m) * y(i,k,m);
        end
    end
end
c_vec = CompDelay(G_workload, CompCostType, CompCostPara);
% Only add cost from nodes with valid computation delay (not NaN or Inf)
valid_nodes = ~isnan(c_vec) & ~isinf(c_vec);
comp_cost = sum(c_vec(valid_nodes) .* G_workload(valid_nodes));

%% Compute Total Utility
total_utility = 0;
for i = 1:N_node
    for k = 1:N_app
        for m = 1:maxN_model
            total_utility = total_utility + Reward(k,m) * (InputRate(i,k) * s(i,k,m));
        end
    end
end

total_request_rate = sum(InputRate(:));

%% Compute Objectives
J_total = -(eta * total_utility - link_cost_total - comp_cost) / total_request_rate;
AvgDelay = (link_cost_total + comp_cost) / total_request_rate;
AvgReward = total_utility / total_request_rate;
%J_no_tun = -(eta * total_utility - link_cost_no_tun - comp_cost) / total_request_rate;
%diffJ = J_no_tun - J_total;
end
