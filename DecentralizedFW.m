function [s_opt, y_opt, Phi_opt, obj_history] = DecentralizedFW(...
    G, InputRate, s_init, y_init, Phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options)
% DecentralizedFW performs a decentralized flow-level Frank–Wolfe optimization
% for joint AI model selection, placement, and routing.
%
%   [s_opt, y_opt, Phi_opt, obj_history] = DecentralizedFW(...
%       InputRate, s_init, y_init, Phi_init, Reward, eta, ...
%       LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
%       Lambda_Node, q_mat, options)
%
%   Inputs:
%       InputRate    - [N_node x N_app] matrix of input rates at nodes.
%       s_init       - [N_node x N_app x maxN_model] initial model selection fractions.
%       y_init       - [N_node x N_app x maxN_model] initial placement fractions.
%       Phi_init     - [N_node x N_app x maxN_model x N_node] initial routing probabilities.
%                      For each node i and each (k,m): sum_j Phi(i,k,m,j) = 1 - y(i,k,m).
%       Reward       - [N_app x maxN_model] matrix of per-request utility values.
%       eta          - Tradeoff weight between delay and utility.
%       LinkCostType - String specifying the link cost function ('linear','queue','taylor').
%       LinkCostPara - [N_node x N_node] matrix of link capacities.
%       CompCostType - String specifying the computation cost function ('linear','queue','taylor').
%       CompCostPara - Vector (length = N_node) of computation capacities.
%       Lambda_Node  - [1 x N_node] vector of node transition rates.
%       q_mat        - [N_node x N_node] movement probability matrix.
%       options      - Structure with fields:
%                        .maxIter      - Maximum number of iterations.
%                        .tol          - Convergence tolerance.
%                        .stepSize_s   - Step size for s updates.
%                        .stepSize_phi - Step size for Phi updates.
%                        .stepSize_y   - Step size for y updates.
%                        .fixPlacement - Boolean; if true, y remains fixed.
%                        .useTunneling - Boolean; if false, tunneling is disabled.
%                        .L_req        - [1 x N_app] vector of request sizes.
%                        .L_res        - [1 x N_app] vector of response sizes.
%                        .L_model      - [N_app x maxN_model] matrix of model hosting costs.
%
%   Outputs:
%       s_opt       - Optimized model selection fractions [N_node x N_app x maxN_model].
%       y_opt       - Optimized placement fractions [N_node x N_app x maxN_model].
%                     (If fixPlacement is true, y_opt equals y_init.)
%       Phi_opt     - Optimized routing probabilities [N_node x N_app x maxN_model x N_node].
%       obj_history - A vector containing the overall objective at each iteration.
%

% Initialization
N_node = size(InputRate,1);
N_app = size(InputRate,2);
maxN_model = size(s_init,3);
global CENTER_NODE;

s = s_init;
y = y_init;
Phi = Phi_init;
obj_history = [0,0];

maxIter = options.maxIter;
tol = options.tol;
stepSize_s = options.stepSize_s;
stepSize_phi = options.stepSize_phi;
if ~options.fixPlacement
    stepSize_y = options.stepSize_y;
end
BlockNodes = computeBlockedNodes(G, CENTER_NODE, false);

for iter = 1:maxIter
    %fprintf('Iteration %d\n', iter);

    %% Step 1: Compute Request and Response Flows
    f_req_cell = CalcRequestFlow(InputRate, s, Phi);   % cell array: [N_app x maxN_model]
    f_res_cell = CalcResponseFlow(f_req_cell);          % response flows as transposes

    G_workload = ComputeNodeWorkload(InputRate, s, f_req_cell, options.L_model);


    %% Step 2: Compute Tunneling Flows & Expected Delay D0
    [D0_cell, f_tun_cell] = CalcTunnelingFlow(InputRate, options.L_req, options.L_res, G_workload,  f_req_cell, f_res_cell, s, y, Phi, Lambda_Node, q_mat, LinkCostType, LinkCostPara, CompCostType, CompCostPara, options.d_AP);

    %% Step 3: Aggregate Flows & Compute Link Delays
    F_req = zeros(N_node);
    F_res = zeros(N_node);
    F_tun = zeros(N_node);
    for k = 1:N_app
        for m = 1:maxN_model
            F_req = F_req + options.L_req(k) * f_req_cell{k,m};
            F_res = F_res + options.L_res(k) * f_res_cell{k,m};
            F_tun = F_tun + options.L_res(k) * f_tun_cell{k,m};
        end
    end

    if ~options.useTunneling
        F_total_alg = F_req + F_res + F_tun;
    else
        F_total_alg = F_req + F_res;
    end


    %% Step 4: Compute Local Data
    % B_ij, m_i^km, etc.


    D_prime = LinkDelayTotalDerivative(F_total_alg, LinkCostType, LinkCostPara);

    if ~options.useTunneling
        d_mat = LinkDelay(F_total_alg, LinkCostType, LinkCostPara);
        d_prime = LinkDelayDerivative(F_total_alg, LinkCostType, LinkCostPara);
        localData = computeLocalData(InputRate, s, D0_cell, Lambda_Node, q_mat, d_mat, F_total_alg, LinkCostType, LinkCostPara, Phi);

        %% Step 5: Compute Accumulated Cost Messages M
        M_cell = computeM(InputRate, s, Phi, localData);
        tau = computeTau(D0_cell, Lambda_Node, q_mat, D_prime);
    else
        localData.r = ones(N_node,N_app,maxN_model);
    end

    %% Step 6: Compute partial J/ partial F_o
    % if consider tunneling, dJdF_o is calcilated by (38)
    if ~options.useTunneling
        dJdF = compute_dJ_dF(D_prime, d_prime, localData.B, Phi, M_cell);
    else
        % if does not consider tunneling, (38) is reduced to previous papers.
        dJdF = D_prime;
    end

    %% Step 7: Compute Exact Delay Sensitivity δ
    delta_cell = computeDelta(Phi, dJdF, y, G_workload, CompCostType, CompCostPara, LinkCostPara, options.L_req, options.L_res);


    %% Step 10: Compute Partial Derivatives
    if ~options.useTunneling
        [grad_s, grad_phi] = computePartialDerivatives(localData, dJdF, delta_cell, tau, Reward, eta, options.L_req,options.L_res);
    else
        [grad_s, grad_phi] = computePartialDerivatives(localData, dJdF, delta_cell, zeros(N_node,N_app,maxN_model), Reward, eta, options.L_req,options.L_res);
    end
    %% Step 11: Frank–Wolfe Update
    if options.fixPlacement % only update s and phi
        [s, Phi] = updateFixedPlacement(BlockNodes, s, Phi, grad_s, grad_phi, stepSize_s, stepSize_phi, LinkCostPara, y);

    elseif ~options.fixPlacement && options.greedyPlacement % LPF + Greedy service placement
        [s, Phi] = updateFixedPlacement(BlockNodes, s, Phi, grad_s, grad_phi, stepSize_s, stepSize_phi, LinkCostPara, y);
        t = computeTotalServedRequests(f_req_cell);  % [N_node x N_app x maxN_model]
        [y, Phi] = updateGreedyPlacementFromT(t, Phi, NodeCap, G,BlockNodes);

    else % DMP-LFW-P
        if iter == 1
            [s, Phi] = updateFixedPlacement(BlockNodes, s, Phi, grad_s, grad_phi, stepSize_s, stepSize_phi, LinkCostPara, y);
            t = computeTotalServedRequests(f_req_cell);  % [N_node x N_app x maxN_model]
            [y, Phi] = updateGreedyPlacementFromT(t, Phi, NodeCap, G,BlockNodes);
        else
            [s, Phi, y] = updateJointPlacement(BlockNodes, s, Phi, y, grad_s, grad_phi, stepSize_s, stepSize_phi, stepSize_y, LinkCostPara, options, NodeCap);

        end
    end
    %% Step 12: Evaluate Overall Objective and Check Convergence
    [J_total, AvgDelay, AvgReward] = CompareObjective(F_req, F_res, F_tun, InputRate, s, y, Reward, eta, CompCostType, CompCostPara, LinkCostType, LinkCostPara);
    %obj_history = [obj_history; J_total];
    obj_history = [obj_history; AvgDelay,AvgReward];
    if iter == 1 || mod(iter,100) == 0
        fprintf('Iteration %d: Objective = %f\n', iter, J_total);
    end
    %% Visualization of flows
    F_total = F_req + F_res + F_tun;
    % if iter == 1
    %     MaxFlow = max(max(F_total));
    % end
    % if iter == 1 || mod(iter,50) == 0
    % 
    %     figure(101); clf;
    %     visualizeNetworkFlow(G, F_total,MaxFlow);
    %     title(sprintf('Network Flow Visualization at Iteration %d', iter));
    %     drawnow;
    % end
    %if iter > 1 && abs(obj_history(end) - obj_history(end-1)) < tol
    %    fprintf('Convergence reached at iteration %d.\n', iter);
    %    break;
    %end

    % check availability
    yphi_sum = ones(N_node,N_app,maxN_model);
    for i = 1:N_node
        for k = 1:N_app
            for m = 1:maxN_model
                yphi_sum(i,k,m) = sum(Phi(i,k,m,:)) + y(i,k,m);
            end
        end
    end
end

s_opt = s;
y_opt = y;
Phi_opt = Phi;
%[s, Phi, y] = updateJointPlacement(BlockNodes, s, Phi, y, grad_s, grad_phi, stepSize_s, stepSize_phi, stepSize_y, LinkCostPara, options, NodeCap);


end
