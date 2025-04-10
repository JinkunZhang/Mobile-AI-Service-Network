function [D0_final, f_tun_final] = CalcTunnelingFlow(InputRate, L_req, L_res, G_workload, f_req_cell, f_res_cell, s, y, Phi, Lambda_Node, q_mat, LinkCostType, LinkCostPara, CompCostType, CompCostPara, d_AP)
% CalcTunnelingFlow computes the expected end-to-end delay D0 (including the effect 
% of tunneling) and the tunneling flows f_tun by iteratively updating them.
%
%   [D0_final, f_tun_final] = CalcTunnelingFlow(InputRate, s, y, Phi, Lambda_Node, ...
%       q_mat, LinkCostType, LinkCostPara, CompCostType, CompCostPara)
%
%   Inputs:
%       InputRate    - [N_node x N_app] matrix of input rates.
%       s            - [N_node x N_app x maxN_model] array of model selection fractions.
%       y            - [N_node x N_app x maxN_model] array of local processing fractions.
%       Phi          - [N_node x N_app x maxN_model x N_node] routing probabilities.
%                      For each node i and each (k,m): sum_j Phi(i,k,m,j) = 1 - y(i,k,m).
%       Lambda_Node  - [1 x N_node] vector of transition rates at each node.
%       q_mat        - [N_node x N_node] matrix where q_mat(i,j) is the probability
%                      that a user at node i moves to node j.
%       LinkCostType - String specifying link cost type ('linear','queue','taylor').
%       LinkCostPara - Link capacities as a vector for valid links.
%                      (We assume here that in main.m, LinkCostPara is defined as a vector.)
%       CompCostType - String specifying computation cost type ('linear','queue','taylor').
%       CompCostPara - Computation capacities at nodes (vector of length N_node).
%
%   Outputs:
%       D0_final    -  [N_app x maxN_model] cell of [N_node x 1] vector, expected end-to-end delays.
%       f_tun_final - [N_app x maxN_model] cell of [N_node x N_node] matrix, tunneling flows.
%
% This function uses an outer iterative procedure with fixed tolerances.

% Set default outer tolerance and maximum iterations.
tol_iter = 1e-5;       
maxIter_iter = 100;

N_node = size(InputRate,1);
N_app = size(InputRate,2);
maxN_model = size(s,3);

%% Convert LinkCostPara vector to full matrix LinkCostMat.
% In main.m, LinkCostPara is defined as a vector for valid edges.
% We assume that invalid links have capacity 0.
% (If your main.m has already done this conversion, you can skip this step.)
% Here we assume a universal definition: LinkCostPara is a full matrix.
LinkCostMat = LinkCostPara;  % For this version, we assume main.m defines LinkCostPara as an [N_node x N_node] matrix.

%% Compute Node Workload and c_vec using CompDelay.
c_vec = CompDelay(G_workload, CompCostType, CompCostPara);

%% Compute Request and Response Flows (independent of tunneling)
%f_req_cell = CalcRequestFlow(InputRate, s, Phi);  % cell array {N_app x maxN_model}
%f_res_cell = CalcResponseFlow(f_req_cell);        % cell array of same size

% Aggregate flows over all (k,m)
F_req = zeros(N_node);
F_res = zeros(N_node);
for k = 1:N_app
    for m = 1:maxN_model
        F_req = F_req + L_req(k) * f_req_cell{k,m};
        F_res = F_res + L_res(k) * f_res_cell{k,m};
    end
end

% Initialize tunneling flows f_tun as zeros.
f_tun = zeros(N_node, N_app, maxN_model, N_node);

%% Outer Iteration: update f_tun and D0 until convergence
D0 = zeros(N_node, N_app, maxN_model);
prev_f_tun_total = zeros(N_node);
for outer_iter = 1:maxIter_iter
    % Aggregate tunneling flows over (k,m) for each link.
    F_tun = zeros(N_node);
    for i = 1:N_node
        for k = 1:N_app
            for m = 1:maxN_model
                F_tun(i,:) = F_tun(i,:) + L_res(k) * squeeze(f_tun(i,k,m,:))';
            end
        end
    end
    % Total flow on each link:
    F_total = F_req + F_res + F_tun;
    
    % Compute link delays using LinkDelay.
    d_mat = LinkDelay(F_total, LinkCostType, LinkCostMat);
    
    % For each (k,m), solve for expected delay D0(:,k,m) using the linear system.
    for k = 1:N_app
        for m = 1:maxN_model
            Phi_km = squeeze(Phi(:, k, m, :));  % [N_node x N_node]
            y_km = squeeze(y(:, k, m));           % [N_node x 1]
            D0_vec = ComputeD0Branch(Phi_km, y_km, c_vec, d_mat, LinkCostMat);
            D0(:,k,m) = D0_vec + 2*d_AP;
        end
    end
    
    % Update tunneling flows f_tun:
    for i = 1:N_node
        for k = 1:N_app
            for m = 1:maxN_model
                r_ikm = InputRate(i,k) * s(i,k,m); % Request rate for (i,k,m)
                for j = 1:N_node
                    if q_mat(i,j) > 0
                        % Tunneling probability as in the paper: p_tun(i,j) = q_mat(i,j) * Lambda_Node(i) * D0(i,k,m)
                        p_tun = q_mat(i,j) * Lambda_Node(i) * D0(i,k,m);
                        % p_tun = q_mat(i,j) * (1 - exp(-Lambda_Node(i) * D0(i,k,m)));
                        f_tun(i,k,m,j) = r_ikm * p_tun;
                    else
                        f_tun(i,k,m,j) = 0;
                    end
                end
            end
        end
    end
    
    % Check convergence: compare aggregated tunneling flows.
    new_f_tun_total = zeros(N_node);
    for i = 1:N_node
        for k = 1:N_app
            for m = 1:maxN_model
                new_f_tun_total(i,:) = new_f_tun_total(i,:) + squeeze(f_tun(i,k,m,:))';
            end
        end
    end
    diff = norm(new_f_tun_total - prev_f_tun_total, 'fro');
    if isnan(diff) || outer_iter == maxIter_iter
        warning('CalcTunnelingFlow did not converge within %d iterations.', maxIter_iter);
        break;
    end

    %fprintf('Outer iteration %d: tunneling flow change = %e\n', outer_iter, diff);
    if diff < tol_iter
        break;
    end
    prev_f_tun_total = new_f_tun_total;
end

%if outer_iter == maxIter_iter
%    warning('CalcTunnelingFlow did not converge within %d iterations.', maxIter_iter);
%end

%% Convert outputs to cell arrays to match format of f_req and f_res.
D0_cell = cell(N_app, maxN_model);
f_tun_cell = cell(N_app, maxN_model);
for k = 1:N_app
    for m = 1:maxN_model
        D0_cell{k,m} = squeeze(D0(:,k,m));            % [N_node x 1] vector.
        f_tun_cell{k,m} = squeeze(f_tun(:,k,m,:));       % [N_node x N_node] matrix.
    end
end

D0_final = D0_cell;
f_tun_final = f_tun_cell;
end
