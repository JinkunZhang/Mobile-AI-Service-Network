function d_prime = LinkDelayDerivative(F, costType, costPara)
% LinkDelayDerivative computes the derivative of the per-link delay with respect
% to the link flow matrix F.
%
%   d_prime = LinkDelayDerivative(F, costType, costPara)
%
%   Inputs:
%       F         - [N_node x N_node] matrix of link flows.
%       costType  - String specifying the cost type. Options:
%                   'linear', 'queue', or 'taylor'
%       costPara  - Link capacity (or cost parameter). Can be a scalar 
%                   or an [N_node x N_node] matrix. For 'queue' and 'taylor',
%                   this represents μ (the service rate).
%
%   Output:
%       d_prime   - [N_node x N_node] matrix of derivatives of the per-link delay
%                   with respect to F.
%
%   Implemented cost functions:
%       1. 'linear'  : If d = costPara .* F, then derivative d' = costPara.
%       2. 'queue'   : If d = 1./(costPara - F) for F < costPara, then
%                      d' = 1./(costPara - F).^2; for F >= costPara, d' = Inf.
%       3. 'taylor'  : If d = 1./costPara + F./(costPara.^2) + (F.^2)./(costPara.^3) + (F.^3)./(costPara.^4),
%                      then derivative d' = 1./(costPara.^2) + 2*F./(costPara.^3) + 3*(F.^2)./(costPara.^4).

switch lower(costType)
    case 'linear'
        % d = costPara .* F   =>   d' = costPara.
        d_prime = costPara;
        
    case 'queue'
        % d = 1./(costPara - F)   =>   d' = 1./(costPara - F).^2.
        d_prime = zeros(size(F));
        for i = 1:size(F,1)
            for j = 1:size(F,2)
                if costPara(i,j) == 0  % Non-link: undefined derivative.
                    d_prime(i,j) = NaN;
                elseif F(i,j) < costPara(i,j)
                    d_prime(i,j) = 1 / ((costPara(i,j) - F(i,j))^2);
                else
                    d_prime(i,j) = Inf;
                    warning('link flow %f exceeds capacity %f on link (%d,%d)', ...
                        F(i,j), costPara(i,j), i, j);
                end
            end
        end
        
    case 'taylor'
        % d = 1./costPara + F./(costPara.^2) + F.^2./(costPara.^3) + F.^3./(costPara.^4)
        % d' = 1./(costPara.^2) + 2*F./(costPara.^3) + 3*(F.^2)./(costPara.^4)
        d_prime = 1./(costPara.^2) + 2*F./(costPara.^3) + 3*(F.^2)./(costPara.^4);
        
    otherwise
        error('Unknown cost type: %s', costType);
end
end
