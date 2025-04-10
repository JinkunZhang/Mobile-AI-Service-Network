function d_mat = LinkDelay(F, costType, costPara)
% LinkDelay computes the per-link delay given the link flow matrix F.
%
%   d_mat = LinkDelay(F, costType, costPara)
%
%   Inputs:
%       F         - [N_node x N_node] matrix of link flows.
%       costType  - String specifying the cost type. Options:
%                   'linear', 'queue', or 'taylor'
%       costPara  - Link capacity (or coefficient) parameter. Can be a scalar 
%                   or an [N_node x N_node] matrix. For 'queue' and 'taylor', 
%                   this represents μ (the service rate).
%
%   Output:
%       d_mat     - [N_node x N_node] matrix of per-link delays.
%
%   Implemented cost functions:
%       1. 'linear'  : d = costPara .* F
%       2. 'queue'   : d = 1./(costPara - F) for F < costPara, else Inf.
%       3. 'taylor'  : 3rd-order Taylor approximation of 1/(costPara - F)

    switch lower(costType)
        case 'linear'
            % Elementwise linear cost.
            d_mat = costPara .* F;
            
        case 'queue'
            % M/M/1 queueing delay: d = 1/(μ - F)
            d_mat = zeros(size(F));
            % Ensure elementwise operations. If F >= costPara, set delay = Inf.
            for i = 1:size(F,1)
                for j = 1:size(F,2)
                    if costPara(i,j) == 0 % if not a link
                        d_mat(i,j) = NaN;
                    elseif F(i,j) < costPara(i,j) % if within capacity
                        d_mat(i,j) = 1 / (costPara(i,j) - F(i,j));
                    else % if outside capacity
                        d_mat(i,j) = Inf;
                        warning('link flow %f exceeds capacity %f on link (%d,%d)', F(i,j), costPara(i,j), i, j);
                    end
                end
            end
            
        case 'taylor'
            % Third-order Taylor expansion of 1/(μ - F) about F = 0.
            % 1/(μ - F) = 1/μ + F/μ^2 + F^2/μ^3 + F^3/μ^4 + ...
            d_mat = 1./costPara + F./(costPara.^2) + (F.^2)./(costPara.^3) + (F.^3)./(costPara.^4);
            
        otherwise
            error('Unknown cost type: %s', costType);
    end
end
