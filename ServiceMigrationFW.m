function [s_opt, y_opt, Phi_opt, obj_history] = ServiceMigrationFW(...
    G, InputRate, s_init, y_init, Phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options)

% ServiceMigrationFW - A variant of decentralized flow-level Frank-Wolfe
% optimization using only in-layer tunneling flow (model-size-based),
% mimicking service migration behavior.

N_node = size(InputRate,1);
N_app = size(InputRate,2);
maxN_model = size(s_init,3);
global CENTER_NODE;

s = s_init;
y = y_init;
Phi = Phi_init;
obj_history = [];
maxIter = options.maxIter;
tol = options.tol;
stepSize_s = options.stepSize_s;
stepSize_phi = options.stepSize_phi;
if ~options.fixPlacement
    stepSize_y = options.stepSize_y;
end

BlockNodes = computeBlockedNodes(G, CENTER_NODE, true);
HopDist = computeHopDistanceToCenter(G, CENTER_NODE);  % for layer info
HopDist = [1;1;1;1;2;2;1;2;3]; % manually set for 9-grid

for iter = 1:maxIter
    %fprintf('Iteration %d\n', iter);

    % Step 1: Compute flows
    f_req_cell = CalcRequestFlow(InputRate, s, Phi);
    f_res_cell = CalcResponseFlow(f_req_cell);
    G_workload = ComputeNodeWorkload(InputRate, s, f_req_cell, options.L_model);

    % Step 2: Compute tunneling flows (model-size, in-layer only)
    [D0_cell, f_tun_cell] = CalcLayeredTunnelingFlow...
    (InputRate, options.L_model, G_workload, ...
    f_req_cell, f_res_cell, s, y, Phi, Lambda_Node, q_mat, LinkCostType, LinkCostPara, ...
    CompCostType, CompCostPara, options.d_AP, HopDist);

    % Step 3: Aggregate flows
    F_req = zeros(N_node);
    F_res = zeros(N_node);
    F_tun = zeros(N_node);
    for k = 1:N_app
        for m = 1:maxN_model
            F_req = F_req + options.L_req(k) * f_req_cell{k,m};
            F_res = F_res + options.L_res(k) * f_res_cell{k,m};
            F_tun = F_tun + options.L_model(k,m) * f_tun_cell{k,m};
        end
    end
    F_total_alg = F_req + F_res;

    % Step 4: Local delay info
    d_mat = LinkDelay(F_total_alg, LinkCostType, LinkCostPara);
    D_prime = LinkDelayTotalDerivative(F_total_alg, LinkCostType, LinkCostPara);
    d_prime = LinkDelayDerivative(F_total_alg, LinkCostType, LinkCostPara);

    localData = computeLocalData(InputRate, s, D0_cell, Lambda_Node, q_mat, d_mat, F_total_alg, LinkCostType, LinkCostPara, Phi);
    M_cell = computeM(InputRate, s, Phi, localData);
    tau = computeTau(D0_cell, Lambda_Node, q_mat, D_prime);

    % Step 5: Derivatives and updates
    %dJdF = compute_dJ_dF(D_prime, d_prime, localData.B, Phi, M_cell);
    dJdF = D_prime;
    delta_cell = computeDelta(Phi, dJdF, y, G_workload, CompCostType, CompCostPara, LinkCostPara, options.L_req, options.L_res);
    %[grad_s, grad_phi] = computePartialDerivatives(localData, dJdF, delta_cell, tau, Reward, eta, options.L_req, options.L_res);
    [grad_s, grad_phi] = computePartialDerivatives(localData, dJdF, delta_cell, zeros(N_node,N_app,maxN_model), Reward, eta, options.L_req,options.L_res);
    

    % Step 6: Frank–Wolfe updates
    if options.fixPlacement
        [s, Phi] = updateFixedPlacement(BlockNodes, s, Phi, grad_s, grad_phi, stepSize_s, stepSize_phi, LinkCostPara, y);
    else
        [s, Phi, y] = updateJointPlacement(BlockNodes, s, Phi, y, grad_s, grad_phi, stepSize_s, stepSize_phi, stepSize_y, LinkCostPara, options, NodeCap);
    end

    % Step 7: Evaluate
    [J_total, ~, ~] = CompareObjective(F_req, F_res, F_tun, InputRate, s, y, Reward, eta, CompCostType, CompCostPara, LinkCostType, LinkCostPara);
    obj_history = [obj_history; J_total];

    if iter == 1 || mod(iter,100) == 0
        fprintf('Objective at iter %d: %f\n', iter, J_total);
    end

    
    % Optional: visualize flows
    % if iter == 1
    %     MaxFlow = max(max(F_total_alg));
    % end
    % if mod(iter, 20) == 0
    %     figure(111); clf;
    %     visualizeNetworkFlow(G, F_total_alg, MaxFlow);
    %     title(sprintf('Service Migration Flow Visualization (Iter %d)', iter));
    %     drawnow;
    % end

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

% Output
s_opt = s;
y_opt = y;
Phi_opt = Phi;
end
