function [s_new, Phi_new, y_new] = updateJointPlacement(BlockNodes, s, Phi, y, grad_s, grad_phi, stepSize_s, stepSize_phi, stepSize_y, LinkCostPara, options, NodeCap)
% updateJointPlacement performs the Frank–Wolfe update for joint optimization,
% updating model selection (s), routing (Phi), and placement (y) variables jointly.
%
%   [s_new, Phi_new, y_new] = updateJointPlacement(BlockNodes, s, Phi, y, grad_s, grad_phi, ...
%         stepSize_s, stepSize_phi, stepSize_y, LinkCostPara, options, NodeCap)
%
%   Inputs:
%       BlockNodes  - Cell array of blocked neighbor indices for each node.
%       s           - [N_node x N_app x maxN_model] current model selection fractions.
%       Phi         - [N_node x N_app x maxN_model x N_node] current routing probabilities.
%       y           - [N_node x N_app x maxN_model] current placement fractions.
%       grad_s      - [N_node x N_app x maxN_model] gradient with respect to s.
%       grad_phi    - [N_node x N_app x maxN_model x N_node] gradient with respect to Phi.
%       stepSize_s  - Scalar step size for s update.
%       stepSize_phi- Scalar step size for Phi update.
%       stepSize_y  - Scalar step size for y update.
%       LinkCostPara- [N_node x N_node] matrix of link capacities.
%       options     - Structure with algorithm options (must include L_model).
%       NodeCap     - [1 x N_node] vector of node capacities.
%
%   Outputs:
%       s_new   - Updated s.
%       Phi_new - Updated Phi.
%       y_new   - Updated placement fractions.
%
global CENTER_NODE;
[N_node, N_app, maxN_model] = size(s);
s_new = s;
Phi_new = Phi;
y_new = y;

%% Step 1: Update s (model selection) using Frank–Wolfe.
for i = 1:N_node
    for k = 1:N_app
        [~, m_star] = min(squeeze(grad_s(i,k,:)));
        descent_s = zeros(maxN_model, 1);
        descent_s(m_star) = 1;
        s_new(i,k,:) = (1 - stepSize_s) * squeeze(s(i,k,:)) + stepSize_s * descent_s;
    end
end

%% Step 2: Joint update for y and Phi.
% For each node i, for each application k, jointly solve for the descent direction d_y by solving an LP,
% then update y and use (1 - d_y) to scale the Phi update.

% Number of decision variables per node for y: numVar = N_app * maxN_model.
numVar = N_app * maxN_model;
L_model_vec = options.L_model(:);  % Column vector of length numVar.

for i = 1:N_node
    % Build objective coefficient vector c for node i: 
    % For each branch (k,m), let c_{i,k,m} = min_{j in allowed} grad_phi(i,k,m,j).
    c = zeros(numVar, 1);
    for k = 1:N_app
        for m = 1:maxN_model
            idx = (k-1)*maxN_model + m;
            valid_j = setdiff(find(LinkCostPara(i,:) > 0), BlockNodes{i});
            if ~isempty(valid_j)
                c(idx) = min(squeeze(grad_phi(i,k,m, valid_j)));
            else
                c(idx) = 0;
            end
        end
    end
    % Solve the LP for node i:
    % maximize   sum_{k,m} c_{i,k,m} * d_y(i,k,m)
    % subject to sum_{k,m} L_model(k,m) * d_y(i,k,m) <= NodeCap(i), and 0 <= d_y <= 1.
    % Because linprog minimizes, we solve:
    % minimize   -c'*x subject to A*x <= NodeCap(i) and 0 <= x <= 1.
    f_obj = -c;
    A = L_model_vec';  % 1 x numVar
    b = NodeCap(i);
    lb = zeros(numVar, 1);
    ub = ones(numVar, 1);
    options_lp = optimoptions('linprog','Display','none');
    %[d_y_sol, ~, exitflag] = linprog(f_obj, A, b, [], [], lb, ub, options_lp);
    
n = length(lb);  % size of your decision variable
try
   [d_y_sol, ~, exitflag] = linprog(f_obj, A, b, [], [], lb, ub, options_lp);
    
    % Additional check in case linprog didn't throw error but returned NaNs
    if any(isnan(d_y_sol))
        warning('linprog returned NaNs. Returning zero solution.');
        d_y_sol = zeros(n,1);
    end

catch ME
    warning('linprog failed: %s. Returning zero solution.');
    d_y_sol = zeros(n,1);
end

    % if exitflag <= 0
    %     warning('LP did not converge for node %d; using uniform descent for y.', i);
    %     d_y_sol = ones(numVar, 1) / numVar;
    % end
    % Reshape d_y_sol to [N_app x maxN_model]
    d_y = reshape(d_y_sol, [N_app, maxN_model]);
    
    % Update y using the LP solution.
    y_new(i, :, :) = (1 - stepSize_y) * squeeze(y(i, :, :)) + stepSize_y * d_y;
    y_new(CENTER_NODE,:,:) = 1;

    % Now, update Phi using the coupling: scale the descent for Phi by (1 - d_y) for each branch.
    for k = 1:N_app
        for m = 1:maxN_model
            valid_j = setdiff(find(LinkCostPara(i,:) > 0), BlockNodes{i});
            if ~isempty(valid_j)
                [~, idx] = min(squeeze(grad_phi(i,k,m, valid_j)));
                j_star = valid_j(idx);
                descent_phi = zeros(N_node,1);
                descent_phi(j_star) = 1;
                Phi_new(i,k,m,:) = (1 - stepSize_phi) * squeeze(Phi(i,k,m,:)) + ...
                    stepSize_phi * (1 - d_y(k, m)) * descent_phi;
            end
        end
    end
end

end
