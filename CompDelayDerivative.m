function comp_deriv = CompDelayDerivative(G_workload, costType, costPara)
% CompDelayDerivative computes the derivative of the per-node computation delay 
% with respect to the node workload.
%
%   comp_deriv = CompDelayDerivative(G_workload, costType, costPara)
%
%   Inputs:
%       G_workload - [N_node x 1] vector of computational workloads.
%       costType   - String specifying the computation cost type: 
%                    'linear', 'queue', or 'taylor'.
%       costPara   - Computation capacity parameter (ν); can be a scalar or vector of length N_node.
%
%   Output:
%       comp_deriv - [N_node x 1] vector of derivatives, i.e. d(c_delay)/d(G_workload).
%
%   See explanation above for details of the derivative computation.

% Ensure costPara is a vector with the same dimensions as G_workload.
costPara = reshape(costPara, size(G_workload));

switch lower(costType)
    case 'linear'
        % For linear cost: c_delay = costPara .* G_workload, derivative = costPara.
        comp_deriv = costPara;
        
    case 'queue'
        % For queueing delay: c_delay = 1/(ν - G_workload), derivative = 1/(ν - G_workload)^2.
        comp_deriv = zeros(size(G_workload));
        for i = 1:length(G_workload)
            if G_workload(i) < costPara(i)
                comp_deriv(i) = 1 / ((costPara(i) - G_workload(i))^2);
            else
                comp_deriv(i) = Inf;
                warning('Workload %f exceeds capacity %f at node %d', G_workload(i), costPara(i), i);
            end
        end
        
    case 'taylor'
        % For the Taylor approximation: 
        % c_delay = 1/ν + G_workload/ν^2 + G_workload^2/ν^3 + G_workload^3/ν^4,
        % so derivative = 1/ν^2 + 2*G_workload/ν^3 + 3*G_workload^2/ν^4.
        comp_deriv = 1./(costPara.^2) + 2*G_workload./(costPara.^3) + 3*(G_workload.^2)./(costPara.^4);
        
    otherwise
        error('Unknown computation cost type: %s', costType);
end
end
