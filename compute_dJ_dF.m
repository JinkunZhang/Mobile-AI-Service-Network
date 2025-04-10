function dJdF = compute_dJ_dF(D_prime, d_prime, B, Phi, M_cell)
% compute_dJ_dF computes the partial derivative of J with respect to the
% outgoing flow F^o, using the updated formula.
%
%   dJdF = compute_dJ_dF(D_prime, d_prime, B, Phi, M_cell)
%
%   Inputs:
%       D_prime  - [N x N] matrix of D'_{ij} (derivative of total delay D_{ij} w.r.t. flow).
%       d_prime  - [N x N] matrix of d'_{ij} (derivative of one-way delay d_{ij} w.r.t. flow).
%       B        - [N x N] matrix of the local derivative aggregate B_{ij}.
%       Phi      - [N x N_app x maxN_model x N] routing probabilities.
%       M_cell   - Cell array of size [N_app x maxN_model], each cell is [N x 1]
%                  representing the accumulated cost M_i^{k,m}.
%
%   Output:
%       dJdF    - [N x N] matrix, where each entry is:
%                 D_prime(i,j) + (d_prime(i,j) / (1 - B(i,j))) * sum_{k,m} Phi(i,k,m,j)*M_i^{k,m}.
%
[N, ~] = size(D_prime);
[N_app, maxN_model] = size(M_cell);
dJdF = zeros(N, N);

for i = 1:N
    for j = 1:N
        % Only consider if the link is valid (i.e., B(i,j) < 1)
        if B(i,j) < 1
            sum_val = 0;
            for k = 1:N_app
                for m = 1:maxN_model
                    % Get Phi(i,k,m,j) and the M value at node i for (k,m)
                    sum_val = sum_val + Phi(i,k,m,j) * M_cell{k,m}(i);
                end
            end
            dJdF(i,j) = D_prime(i,j) + (d_prime(i,j) / (1 - B(i,j))) * sum_val;
        else
            dJdF(i,j) = Inf; % or NaN, if the denominator is zero
            warning('B exceeding 1 on link (%d,%d)',i,j);
        end
    end
end

end
