function q_mat = computeMobilityMatrix(G, initOption, visualize)
% computeMobilityMatrix returns the mobility probability matrix q_mat.
%
%   q_mat = computeMobilityMatrix(G, initOption, visualize)
%
%   Inputs:
%       G         - MATLAB graph object. Assumes G.Nodes contains fields 'X' and 'Y'
%                   (if not, sorting by node index will be used).
%       initOption- String: either 'uniform' or 'cyclic'.
%       visualize - (Optional) Boolean flag; if true, the function visualizes
%                   the mobility flows. Default is false.
%
%   Output:
%       q_mat     - [N_node x N_node] matrix where q_mat(i,j) is the probability 
%                   that a user at node i moves to node j.
%
%   For 'uniform': For each node i, assign equal probability to all outgoing neighbors.
%
%   For 'cyclic': First sort nodes by their (X,Y) coordinates (or by index if not available).
%       Then, for each node i (in sorted order), choose as the next hop the node that
%       is next in the sorted order and is a neighbor of i. If not found, choose the
%       neighbor with the smallest index.
%
%   Additionally, if visualize is true, the function overlays arrows on the graph
%   (with thickness proportional to the probability).
%
if nargin < 3 || isempty(visualize)
    visualize = false;
end
global CENTER_NODE;
N_node = numnodes(G);
q_mat = zeros(N_node);

switch lower(initOption)
    case 'uniform'
        % For each node, get its neighbors and assign equal probability.
        for i = 1:N_node
            nbrs = neighbors(G, i);
            if ~isempty(nbrs)
                prob = 1 / length(nbrs);
                q_mat(i, nbrs) = prob;
            end
        end
        
    case 'cyclic'
        % Try to sort nodes using coordinates if available.
        if isfield(G.Nodes, 'X') && isfield(G.Nodes, 'Y')
            coords = [G.Nodes.X, G.Nodes.Y];
            [~, sortIdx] = sortrows(coords);
        else
            sortIdx = (1:N_node)';
        end
        
        % For each node in sorted order, assign its next-hop.
        for ii = 1:N_node
            i = sortIdx(ii);
            nbrs = neighbors(G, i);
            if isempty(nbrs)
                continue;
            end
            if ii < N_node
                ideal_next = sortIdx(ii+1);
            else
                ideal_next = sortIdx(1);
            end
            
            if any(nbrs == ideal_next)
                nextHop = ideal_next;
            else
                nextHop = min(nbrs);
            end
            
            q_mat(i, nextHop) = 1;
        end

       case 'outward'
        % Determine center node (assume min eccentricity)
        D = distances(G);
        dist_to_center = D(:, CENTER_NODE);

        for i = 1:N_node
            nbrs = neighbors(G, i);
            if isempty(nbrs)
                continue;
            end
            % Outward neighbors = neighbors farther from center
            outward_nbrs = nbrs(dist_to_center(nbrs) > dist_to_center(i));
            if ~isempty(outward_nbrs)
                prob = 1 / length(outward_nbrs);
                q_mat(i, outward_nbrs) = prob;
            else
                % fallback: cyclic rule (move to min-indexed neighbor)
                q_mat(i, min(nbrs)) = 1;
            end
        end
        
    otherwise
        error('Unknown mobility option: %s. Use ''uniform'' or ''cyclic''.', initOption);
end

if visualize
    figure;
    % Plot the original graph in light gray.
    h = plot(G, 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1, 'NodeLabel', {});
    hold on;
    
    % Retrieve node coordinates.
    if isfield(G.Nodes, 'X') && isfield(G.Nodes, 'Y')
        X = G.Nodes.X;
        Y = G.Nodes.Y;
    else
        p = plot(G);
        X = p.XData;
        Y = p.YData;
        delete(p);
    end
    
    % Overlay arrows for each nonzero mobility probability.
    for i = 1:N_node
        for j = 1:N_node
            if q_mat(i,j) > 0
                % Set arrow thickness proportional to q_mat(i,j).
                lw = 2 * q_mat(i,j);
                dx = X(j) - X(i);
                dy = Y(j) - Y(i);
                quiver(X(i), Y(i), dx, dy, 0, 'LineWidth', lw, 'Color', 'r', 'MaxHeadSize', 0.5);
            end
        end
    end
    title(sprintf('Mobility Matrix Visualization (%s mode)', initOption));
    hold off;
end
end
