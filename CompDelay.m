function c_delay = CompDelay(G_workload, costType, costPara)
% CompDelay computes the per-node computation delay given the node workload.
%
%   c_delay = CompDelay(G_workload, costType, costPara)
%
%   Inputs:
%       G_workload - [N_node x 1] vector of computational workloads.
%       costType   - String specifying the cost type. Options:
%                    'linear', 'queue', or 'taylor'
%       costPara   - Computation capacity parameter (ν) or coefficient.
%                    Can be a scalar or vector of length N_node.
%
%   Output:
%       c_delay    - [N_node x 1] vector of per-node computation delays.
%
%   Implemented cost functions:
%       1. 'linear'  : c_delay = costPara .* G_workload
%       2. 'queue'   : c_delay = 1./(costPara - G_workload), for G_workload < costPara, else Inf.
%       3. 'taylor'  : 3rd-order Taylor approximation of 1/(costPara - G_workload).
    costPara = reshape(costPara,size(G_workload));
    switch lower(costType)
        case 'linear'
            % Linear cost: elementwise multiplication.
            c_delay = costPara .* G_workload;
            
        case 'queue'
            % Queueing delay: c_delay = 1/(ν - G_workload)
            c_delay = zeros(size(G_workload));
            for i = 1:length(G_workload)
                if G_workload(i) < costPara(i)
                    c_delay(i) = 1 / (costPara(i) - G_workload(i));
                else
                    c_delay(i) = Inf;
                    warning('workload %f exceeds capacity %f on node %d', G_workload(i), costPara(i), i);
                end
            end
            
        case 'taylor'
            % Third-order Taylor approximation of 1/(ν - G_workload)
            c_delay = 1./costPara + G_workload./(costPara.^2) + (G_workload.^2)./(costPara.^3) + (G_workload.^3)./(costPara.^4);
            
        otherwise
            error('Unknown cost type: %s', costType);
    end
end
