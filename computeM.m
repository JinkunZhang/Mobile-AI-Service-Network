function M_cell = computeM(InputRate, s, Phi, localData)
% computeM computes the accumulated cost M for each node, application, and model.
%
%   M_cell = computeM(InputRate, s, y, Phi, localData)
%
%   Inputs:
%       InputRate   - [N_node x N_app] matrix of input rates.
%       s           - [N_node x N_app x maxN_model] array of model selection fractions.
%       y           - [N_node x N_app x maxN_model] array of local processing fractions.
%       Phi         - [N_node x N_app x maxN_model x N_node] routing probabilities.
%                     Here, Phi(i,k,m,j) is the fraction of flow sent from node i to node j.
%                     In the paper, the recursion is defined over incoming flows.
%       localData   - Structure output from computeLocalData, containing local cost m.
%                     localData.m is [N_node x N_app x maxN_model].
%
%   Output:
%       M_cell      - Cell array of size [N_app x maxN_model], where each cell is an
%                     [N_node x 1] vector representing the accumulated cost M for that (k,m).
%
%   Method:
%     For each application k and model m, we solve:
%         (I - (Phi(:,:,k,m))^T) * M = m,
%     where m = localData.m(:,k,m).
%
%   This is equivalent to the recursive messaging protocol described in the paper.
%

N_node = size(InputRate,1);
N_app = size(InputRate,2);
maxN_model = size(s,3);

M_cell = cell(N_app, maxN_model);
I = eye(N_node);

for k = 1:N_app
    for m_idx = 1:maxN_model
        % Extract the routing matrix for application k, model m.
        % In our data, Phi(i,k,m,j) is the fraction from i to j.
        % The recursion in the paper is defined as:
        %   M(i) = m(i) + sum_l Phi(l,k,m,i)*M(l),
        % which is equivalent to: (I - Phi(:,:,k,m)' ) * M = m.
        Phi_km = squeeze(Phi(:, k, m_idx, :));   % [N_node x N_node]
        m_vec = squeeze(localData.m(:, k, m_idx)); % [N_node x 1]
        M_cell{k, m_idx} = (I - Phi_km') \ m_vec;
    end
end

end
