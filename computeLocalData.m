function localData = computeLocalData(InputRate, s, D0_cell, Lambda_Node, q_mat, d_mat, F_total, LinkCostType, LinkCostPara, Phi)
% computeLocalData computes local quantities for derivative calculations.
%
%   localData = computeLocalData(InputRate, s, D0_cell, Lambda_Node, q_mat, d_mat, F_total, ...
%                                 LinkCostType, LinkCostPara, Phi)
%
%   Inputs:
%       InputRate   - [N_node x N_app] matrix of input rates.
%       s           - [N_node x N_app x maxN_model] array of model selection fractions.
%       D0_cell     - Cell array of size [N_app x maxN_model], each cell is [N_node x 1]
%                     representing the expected end-to-end delays for that (k,m).
%       Lambda_Node - [1 x N_node] vector of node transition rates.
%       q_mat       - [N_node x N_node] matrix of movement probabilities.
%       d_mat       - [N_node x N_node] matrix of one-way link delays.
%       F_total     - [N_node x N_node] matrix of total flows on each link.
%       LinkCostType- String specifying the link cost type ('linear','queue','taylor').
%       LinkCostPara- [N_node x N_node] matrix of link capacities (nonzero only for valid links).
%       Phi         - [N_node x N_app x maxN_model x N_node] routing probabilities.
%
%   Output:
%       localData - Structure with fields:
%                      r          - [N_node x N_app x maxN_model] effective request rates.
%                      exp_factor - [N_node x N_app x maxN_model] exp(-Lambda_Node(i)*D0(i,k,m)).
%                      out_cost   - [N_node x 1] vector, where for each node i:
%                                   out_cost(i) = sum_{j in nbr(i)} d_mat(i,j)*q_mat(i,j), 
%                                   summing only over valid links (LinkCostPara(i,j)>0).
%                      m          - [N_node x N_app x maxN_model] local cost contributions,
%                                   m(i,k,m) = Lambda_Node(i)*r(i,k,m)*exp_factor(i,k,m)*out_cost(i).
%                      B          - [N_node x N_node] matrix, where:
%                                   B(i,j) = Lambda_Node(i)*q_mat(i,j)*d_prime(i,j)* ...
%                                             sum_{k=1}^{N_app} sum_{m=1}^{maxN_model}[ r(i,k,m)*Phi(i,k,m,j)*exp_factor(i,k,m) ].
%                      D0         - The input expected delay cell array.
%
%   This function calls LinkDelayDerivative to compute d_prime.
%
N_node = size(InputRate,1);
N_app = size(InputRate,2);
maxN_model = size(s,3);

%% Compute effective request rate and exp_factor.
r_eff = zeros(N_node, N_app, maxN_model);
exp_factor = zeros(N_node, N_app, maxN_model);
for i = 1:N_node
    for k = 1:N_app
        for m = 1:maxN_model
            r_eff(i,k,m) = InputRate(i,k) * s(i,k,m);
            % Retrieve expected delay for node i, app k, model m from D0_cell.
            exp_factor(i,k,m) = exp(-Lambda_Node(i) * D0_cell{k,m}(i));
        end
    end
end

%% Compute out_cost for each node (only over valid links based on LinkCostPara).
out_cost = zeros(N_node, 1);
for i = 1:N_node
    % Identify neighbors j where LinkCostPara(i,j) > 0.
    validNbrs = find(LinkCostPara(i,:) > 0);
    out_cost(i) = sum(d_mat(i, validNbrs) .* q_mat(i, validNbrs));
end

%% Compute local cost contribution m.
m_val = zeros(N_node, N_app, maxN_model);
for i = 1:N_node
    for k = 1:N_app
        for m_idx = 1:maxN_model
            m_val(i,k,m_idx) = Lambda_Node(i) * r_eff(i,k,m_idx) * exp_factor(i,k,m_idx) * out_cost(i);
        end
    end
end

%% Compute the derivative d_prime using LinkDelayDerivative.
d_prime = LinkDelayDerivative(F_total, LinkCostType, LinkCostPara);

%% Compute B matrix.
B = zeros(N_node);
for i = 1:N_node
    for j = 1:N_node
        % Only consider valid links based on LinkCostPara.
        if LinkCostPara(i,j) > 0
            sum_val = 0;
            for k = 1:N_app
                for m_idx = 1:maxN_model
                    sum_val = sum_val + r_eff(i,k,m_idx) * Phi(i,k,m_idx,j) * exp_factor(i,k,m_idx);
                end
            end
            B(i,j) = Lambda_Node(i) * q_mat(i,j) * d_prime(i,j) * sum_val;
        else
            B(i,j) = 0;
        end
    end
end

%% Package computed values into the output structure.
localData.r = r_eff;
localData.exp_factor = exp_factor;
localData.out_cost = out_cost;
localData.m = m_val;
localData.B = B;
localData.D0 = D0_cell;
end
