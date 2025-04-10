function f_res = CalcResponseFlow(f_req)
% CalcResponseFlow computes the response flow on each link for each (k,m).
%
%   f_res = CalcResponseFlow(f_req)
%
%   Input:
%       f_req : cell array of request flows {N_app x maxN_model}, each cell is
%               an [N_node x N_node] matrix.
%
%   Output:
%       f_res : cell array of response flows (same size as f_req).
%
% We assume that the response flow retraces the request path in reverse.
for k = 1:size(f_req,1)
    for m = 1:size(f_req,2)
        f_res{k, m} = f_req{k, m}';
    end
end
end
