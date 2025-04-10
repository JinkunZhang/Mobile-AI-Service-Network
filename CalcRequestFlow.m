function f_req = CalcRequestFlow(InputRate, s, Phi)
% CalcRequestFlow computes the request flow on each link for each (k,m).
%
%   f_req = CalcRequestFlow(InputRate, s, Phi)
%
%   Inputs:
%       InputRate : [N_node x N_app] matrix with input rates at each node for each app.
%       s         : [N_node x N_app x maxN_model] array with model selection fractions.
%       Phi       : [N_node x N_app x maxN_model x N_node] array with routing fractions.
%       L_req     : [1 x N_app] vector, request size for each application.
%
%   Output:
%       f_req : a cell array of size {N_app x maxN_model}, where each cell is an
%               [N_node x N_node] matrix representing the flow on each link.
%
% For each application k and model m, we compute the effective input rate:
%   r_{i,k,m} = InputRate(i,k) * s(i,k,m)
% and then compute the link flow using:
%   t_vec = (I - Phi(:,:,k,m)') \ r  and  f = diag(t_vec) * Phi(:,:,k,m)

N_node = size(InputRate, 1);
N_app = size(InputRate, 2);
maxN_model = size(s, 3);

f_req = cell(N_app, maxN_model);

for k = 1:N_app
    for m = 1:maxN_model
        % r_km is a column vector [N_node x 1]
        r_km = InputRate(:, k) .* squeeze(s(:, k, m));
        % Extract the Phi matrix for app k and model m: [N_node x N_node]
        Phi_km = squeeze(Phi(:, k, m, :));
        I = eye(N_node);
        % Compute t_vec = (I - Phi_km') \ r_km.
        if rcond(I - Phi_km') < 1e-12
            warning('Matrix (I - Phi(:,:,%d,%d)'' ) is nearly singular.', k, m);
        end
        t_vec = (I - Phi_km') \ r_km;
        f_mat = diag(t_vec) * Phi_km;
        f_req{k, m} = f_mat;
    end
end
end
