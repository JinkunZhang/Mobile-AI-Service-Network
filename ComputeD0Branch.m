function D0_vec = ComputeD0Branch(Phi_km, y_km, c_vec, d_mat, LinkCostMat)
% ComputeD0Branch computes the expected end-to-end delay D0 for a given branch (k,m).
%
%   D0_vec = ComputeD0Branch(Phi_km, y_km, c_vec, d_mat, LinkCostMat)
%
%   Inputs:
%       Phi_km      - [N_node x N_node] routing probability matrix for a given (k,m),
%                     where Phi_km(i,j) is the fraction of flow forwarded from node i to node j.
%       y_km        - [N_node x 1] vector of local processing fractions for branch (k,m).
%       c_vec       - [N_node x 1] vector of computation delays at each node.
%       d_mat       - [N_node x N_node] matrix of one-way link delays.
%       LinkCostMat - [N_node x N_node] matrix of link capacities (nonzero for valid links).
%
%   Output:
%       D0_vec      - [N_node x 1] vector of expected end-to-end delays for the branch.
%
%   The function implements the following steps:
%       1. For each node i, compute:
%             b(i) = sum_{j in valid neighbors} Phi_km(i,j) * ( d_mat(i,j) + d_mat(j,i) )
%          where "valid neighbors" are j such that LinkCostMat(i,j) > 0.
%       2. Form the right-hand side vector:
%             r_vec = y_km .* c_vec + b.
%       3. Solve the linear system:
%             (I - Phi_km) * D0_vec = r_vec.
%

N_node = size(Phi_km, 1);
b = zeros(N_node, 1);

for i = 1:N_node
    % Find valid neighbors: those j with nonzero link capacity.
    validNbrs = find(LinkCostMat(i,:) > 0);
    % Sum the contributions for node i.
    b(i) = sum( Phi_km(i, validNbrs) .* ( d_mat(i, validNbrs) + d_mat(validNbrs, i)' ) );
end

% Form the right-hand side.
r_vec = y_km .* c_vec + b;
I = eye(N_node);

% Solve for D0.
D0_vec = (I - Phi_km) \ r_vec;
end
