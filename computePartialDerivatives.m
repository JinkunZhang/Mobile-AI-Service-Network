function [grad_s, grad_phi] = computePartialDerivatives(localData, dJdF, delta_cell, tau, Reward, eta, L_req,L_res)
% computePartialDerivativesExact computes the partial derivatives of the overall
% objective with respect to the decision variables s and Phi.
%
%   [grad_s, grad_phi] = computePartialDerivativesExact(localData, M_cell, delta_cell, tau, Reward, eta, D_prime, d_prime, Phi)
%
%   Inputs:
%       localData  - Structure with fields:
%                    .r: [N_node x N_app x maxN_model] effective request rates.
%                    .B: [N_node x N_node] matrix from local data.
%       M_cell     - Cell array of size [N_app x maxN_model], each cell is [N_node x 1]
%                    representing accumulated cost M_i^{k,m}.
%       delta_cell - Cell array of size [N_app x maxN_model], each cell is [N_node x 1]
%                    representing the exact delay sensitivity δ for (k,m).
%       tau        - [N_node x N_app x maxN_model] array of tunneling delay contributions.
%       Reward     - [N_app x maxN_model] matrix of per-request utility values.
%       eta        - Scalar tradeoff weight.
%       D_prime    - [N_node x N_node] matrix of D'_{ij} (total delay derivative).
%       d_prime    - [N_node x N_node] matrix of d'_{ij} (one-way delay derivative).
%       Phi        - [N_node x N_app x maxN_model x N_node] routing probabilities.
%
%   Outputs:
%       grad_s   - [N_node x N_app x maxN_model] gradient with respect to s.
%       grad_phi - [N_node x N_app x maxN_model x N_node] gradient with respect to Phi.
%
%   The computations follow the formulas:
%
%       (1) Compute dJdF, the partial derivative of J with respect to the outgoing flow:
%           dJdF(i,j) = D_prime(i,j) + (d_prime(i,j)/(1 - B(i,j))) * sum_{k,m} [Phi(i,k,m,j)*M_i^{k,m}].
%
%       (2) Then, for each node i, application k, model m:
%           grad_s(i,k,m) = localData.r(i,k,m) * ( delta_cell{k,m}(i) + tau(i,k,m) - eta*Reward(k,m) ),
%           grad_phi(i,k,m,j) = localData.r(i,k,m) * ( dJdF(i,j) + delta_cell{k,m}(j) ).
%
[N_node, N_app, maxN_model] = size(localData.r);

grad_s = zeros(N_node, N_app, maxN_model);
grad_phi = zeros(N_node, N_app, maxN_model, N_node);

% For each (k,m), update gradients.
for k = 1:N_app
    for m = 1:maxN_model
        delta_vec = delta_cell{k,m}; % [N_node x 1]
        for i = 1:N_node
            % Use effective request rate r(i,k,m) from localData.
            r_val = localData.r(i,k,m);
            % Compute grad_s:
                grad_s(i,k,m) = r_val * ( delta_vec(i) + tau(i,k,m) - eta * Reward(k,m) );
   
            
            % Compute grad_phi for each outgoing neighbor j.
            for j = 1:N_node
                grad_phi(i,k,m,j) = r_val * ( L_req(k)*dJdF(i,j) + L_res(k)*dJdF(j,i) + delta_vec(j) );
            end
        end
    end
end

end
