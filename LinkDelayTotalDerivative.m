function D_prime = LinkDelayTotalDerivative(F, costType, costPara)
% LinkDelayTotalDerivative computes the derivative of the total link delay D with respect to F.
%
%   D_prime = LinkDelayTotalDerivative(F, costType, costPara)
%
%   Inputs:
%       F         - [N x N] matrix of link flows.
%       costType  - String specifying the cost type. Options:
%                   'linear', 'queue', or 'taylor'
%       costPara  - [N x N] matrix of link capacities (or coefficients).
%
%   Output:
%       D_prime   - [N x N] matrix where each entry is:
%                   D' = d(F) + F .* d'(F)
%
%   For the different cost models, we have:
%
%   1. Linear:
%         d = costPara .* F   (elementwise multiplication)
%         D = F .* d = costPara .* F.^2
%         D' = 2 * costPara .* F
%
%   2. Queue:
%         d = 1./(costPara - F)      for F < costPara, else Inf
%         D = F ./ (costPara - F)
%         D' = costPara ./ (costPara - F).^2
%
%   3. Taylor (3rd-order approximation):
%         d = 1./costPara + F./(costPara.^2) + F.^2./(costPara.^3) + F.^3./(costPara.^4)
%         D = F .* d = F./costPara + F.^2./(costPara.^2) + F.^3./(costPara.^3) + F.^4./(costPara.^4)
%         D' = 1./costPara + 2*F./(costPara.^2) + 3*F.^2./(costPara.^3) + 4*F.^3./(costPara.^4)
%

switch lower(costType)
    case 'linear'
        % Compute d and d'
        d_mat = costPara .* F;
        d_prime = costPara; % derivative of d with respect to F is costPara (constant)
        D_prime = d_mat + F .* d_prime;  % = costPara .* F + costPara .* F = 2*costPara.*F
        
    case 'queue'
        % For queue: d = 1./(costPara - F) for F < costPara, else Inf.
        % d' = 1./(costPara - F).^2.
        d_mat = zeros(size(F));
        d_prime = zeros(size(F));
        D_prime = zeros(size(F));
        [N, ~] = size(F);
        for i = 1:N
            for j = 1:N
                if costPara(i,j) == 0
                    d_mat(i,j) = NaN;
                    d_prime(i,j) = NaN;
                    D_prime(i,j) = NaN;
                elseif F(i,j) < costPara(i,j)
                    d_mat(i,j) = 1 / (costPara(i,j) - F(i,j));
                    d_prime(i,j) = 1 / ((costPara(i,j) - F(i,j))^2);
                    D_prime(i,j) = d_mat(i,j) + F(i,j)*d_prime(i,j);
                else
                    d_mat(i,j) = Inf;
                    d_prime(i,j) = Inf;
                    D_prime(i,j) = Inf;
                    warning('Flow %f exceeds capacity %f on link (%d,%d)', F(i,j), costPara(i,j), i, j);
                end
            end
        end
        
    case 'taylor'
        % Third-order Taylor expansion.
        % d = 1./costPara + F./(costPara.^2) + F.^2./(costPara.^3) + F.^3./(costPara.^4)
        % d' = 1./(costPara.^2) + 2*F./(costPara.^3) + 3*F.^2./(costPara.^4)
        d_mat = 1./costPara + F./(costPara.^2) + (F.^2)./(costPara.^3) + (F.^3)./(costPara.^4);
        d_prime = 1./(costPara.^2) + 2*F./(costPara.^3) + 3*(F.^2)./(costPara.^4);
        D_prime = d_mat + F .* d_prime;
        
    otherwise
        error('Unknown cost type: %s', costType);
end

end
