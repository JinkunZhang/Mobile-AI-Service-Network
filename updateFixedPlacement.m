function [s_new, Phi_new] = updateFixedPlacement(BlockNodes, s, Phi, grad_s, grad_phi, stepSize_s, stepSize_phi, LinkCostPara, y)
% updateFixedPlacement performs the Frank–Wolfe update for fixed placement.
% Only the model selection (s) and routing (Phi) variables are updated;
% the placement variables y remain unchanged.
%
%   [s_new, Phi_new] = updateFixedPlacement(s, Phi, grad_s, grad_phi, stepSize_s, stepSize_phi, LinkCostPara, y)
%
%   Inputs:
%       s           - [N_node x N_app x maxN_model] current model selection fractions.
%       Phi         - [N_node x N_app x maxN_model x N_node] current routing probabilities.
%       grad_s      - [N_node x N_app x maxN_model] gradient with respect to s.
%       grad_phi    - [N_node x N_app x maxN_model x N_node] gradient with respect to Phi.
%       stepSize_s  - Scalar step size for s update.
%       stepSize_phi- Scalar step size for Phi update.
%       LinkCostPara- [N_node x N_node] matrix of link capacities; used here to determine allowed links.
%       y           - [N_node x N_app x maxN_model] placement fractions (remain fixed).
%
%   Outputs:
%       s_new       - Updated s.
%       Phi_new     - Updated Phi.
%
[N_node, N_app, maxN_model] = size(s);
s_new = s;
Phi_new = Phi;

% Update s: for each node and each application, select the model with minimal gradient.
for i = 1:N_node
    for k = 1:N_app
        % Find m* that minimizes grad_s(i,k,:)
        [~, m_star] = min(squeeze(grad_s(i,k,:)));
        % Create a descent direction: a unit vector in the m_star coordinate.
        descent_s = zeros(maxN_model, 1);
        descent_s(m_star) = 1;
        % Frank–Wolfe convex combination update.
        s_new(i,k,:) = (1 - stepSize_s) * squeeze(s(i,k,:)) + stepSize_s * descent_s;
    end
end

% Update Phi: for each node, app, and model, choose the outgoing neighbor with minimal gradient.
for i = 1:N_node
    for k = 1:N_app
        for m = 1:maxN_model
            % Allowed outgoing neighbors: where LinkCostPara(i,j) > 0.
           valid_j = setdiff(find(LinkCostPara(i,:) > 0), BlockNodes{i});
            if ~isempty(valid_j)
                [~, idx] = min(squeeze(grad_phi(i,k,m, valid_j)));
                j_star = valid_j(idx);
                descent_phi = zeros(N_node, 1);
                descent_phi(j_star) = 1;
                % For nodes that host the model (y nearly 1), we leave Phi at zero.
                if abs(y(i,k,m) - 1) < 1e-3
                    Phi_new(i,k,m,:) = 0;
                else
                    Phi_new(i,k,m,:) = (1 - stepSize_phi) * squeeze(Phi(i,k,m,:)) + stepSize_phi * descent_phi;
                end
            end
        end
    end
end
end
