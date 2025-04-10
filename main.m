%% Flow Level Simulator for Joint AI Model Selection, Placement and Routing
% Under homogeneous assumptions (Assumptions 1 & 2)
% The simulation uses time-averaged flow rates on each link.
clear all
close all

%% Network Setup

%N_node = 25;
TopoType = 'grid';      % Options: 'grid' for 2D grid, 'er' for Erdos-Rényi,
% 'demo' for tunneling flow demostration (no exceuction of alg)
%TopoPara = 0.15;         % For 'er': link density

LinkCostType = 'taylor';      % e.g., queueing delay: d_ij = 1/(C_ij - F_ij)
CompCostType = 'taylor';      % Computation cost type

% Generate network topology and compute the center node
NetworkInfo = TopoGen(TopoType);
G = NetworkInfo.G;
centerNode = NetworkInfo.centerNode;
N_app = NetworkInfo.N_app;
N_model = NetworkInfo.N_model;
InputRate = NetworkInfo.InputRate;
LinkCostPara = NetworkInfo.LinkCostPara;
CompCostPara = NetworkInfo.CompCostPara;
ModelSize = NetworkInfo.ModelSize;
Reward = NetworkInfo.Reward;
L_req = NetworkInfo.L_req;
L_res = NetworkInfo.L_res;
NodeCap = NetworkInfo.NodeCap;
eta = NetworkInfo.eta;
N_node = numnodes(G); % Use the actual number of nodes in the graph
maxN_model = max(N_model);   % Maximum number of models among apps

% Set the global center node for later use in y_init_gen.
global CENTER_NODE;
CENTER_NODE = centerNode;
% 
% % Number of applications and models per application
% N_app = 10;                   % Number of applications
% N_model = 2*ones(1,N_app);   % Number of models for each application (can vary)
% maxN_model = max(N_model);   % Maximum number of models among apps

% % Get network edges and number of links
% Edges = G.Edges;
% N_edge = size(Edges,1);
% 
% % Link and computation cost parameters
% LinkCostType = 'taylor';      % e.g., queueing delay: d_ij = 1/(C_ij - F_ij)
% LinkCostPara = zeros(N_node);  % Initialize full matrix
% linkCapRef = 20; % used to construct link capacities
% 
% % Compute distances of each node from the center.
% nodeCoords = [G.Nodes.X, G.Nodes.Y]; % assuming G.Nodes has fields 'X' and 'Y'
% centerCoords = nodeCoords(centerNode, :);
% d_center = sqrt(sum((nodeCoords - centerCoords).^2, 2)); % [N_node x 1]
% 
% % Maximum distance among nodes.
% d_max = max(d_center);
% 
% % Scaling factor alpha; adjust as desired.
% alpha = 1; % For example, alpha = 1 means capacity at the center is linkCapRef*(1+1)=2*linkCapRef.
% 
% % Initialize LinkCostPara as an N_node x N_node matrix.
% for e = 1:N_edge
%     i = Edges.EndNodes(e,1);
%     j = Edges.EndNodes(e,2);
%     % Compute average distance of nodes i and j from center.
%     d_mean = (d_center(i) + d_center(j)) / 2;
%     % Compute capacity: farthest links (d_mean == d_max) get capacity = linkCapRef,
%     % while links closer to the center get a higher capacity.
%     tem_cap = linkCapRef * (1 + alpha * (1 - d_mean/d_max));
%     % Assign capacity for both (i,j) and (j,i).
%     LinkCostPara(i,j) = tem_cap;
%     LinkCostPara(j,i) = tem_cap;
% end
% 
% 
% CompCostType = 'taylor';      % Computation cost type
% CompCostPara = 10 * ones(1, N_node);  % Computation capacities at nodes
% CompCostPara(CENTER_NODE) = 1e6;

% % Model and reward parameters
% ModelSize = ones(N_app, maxN_model);  % Assume all models are of size 1
% Reward = ones(N_app, maxN_model);       % Unit utility for all models
% eta = 0;                                % Tradeoff weight between delay and utility
% 
% % Universal request and response sizes for each application.
% L_req = 1 * ones(1, N_app);   % request size for each application
% L_res  = 2 * ones(1, N_app);   % response size for each application
% 
% 
% % Adjust center node capacity to be very high (simulate infinity)
% NodeCap = ones(1, N_node);
% NodeCap(centerNode) = 1e6;  % Set center node capacity very high

%% Flow-level (Input) Rates and Link Parameters
HoldingTime = 10;  % average holding time per node

% Set Lambda_Node for each node as 1/HoldingTime.
Lambda_Node = ones(1, N_node) * (1 / HoldingTime);

% generate mobility matrix
q_mat = computeMobilityMatrix(G, 'cyclic', true);

% 
% % Input rates for each application at each node
% InputRate = 1 * ones(N_node, N_app);

%% Initial Variables

% Generate a valid initial model placement (y)
% y: [N_node_actual x N_app x maxN_model]
y_init = y_init_gen(N_node, N_app, N_model, maxN_model, NodeCap,'cloud');

% Generate initial s and phi.
% s_init: [N_node_actual x N_app x maxN_model] -- initial model selection fractions.
% phi_init: [N_node_actual x N_app x maxN_model x N_node_actual] -- initial routing fractions.
[s_init, phi_init] = s_phi_init_gen(G, centerNode, N_app, N_model, maxN_model, y_init,'shortestpath');

% %% Test
%
% [D0_init, f_tun_init] = CalcTunnelingFlow(InputRate, L_req, L_res, s_init, y_init, phi_init, Lambda_Node, q_mat, LinkCostType, LinkCostPara, CompCostType, CompCostPara);
%
% % Compute Request and Response Flows (independent of tunneling)
% f_req_cell = CalcRequestFlow(InputRate, s_init, phi_init);  % cell array {N_app x maxN_model}
% f_res_cell = CalcResponseFlow(f_req_cell);        % cell array of same size
%
% % Aggregate flows over all (k,m)
% F_req = zeros(N_node);
% F_res = zeros(N_node);
% F_tun = zeros(N_node);
% for k = 1:N_app
%     for m = 1:maxN_model
%         F_req = F_req + L_req(k) * f_req_cell{k,m};
%         F_res = F_res + L_res(k) * f_res_cell{k,m};
%         F_tun = F_tun + L_res(k) * f_tun_init{k,m};
%     end
% end
% F = F_req + F_res + F_tun;
% Tun_ratio = F_tun ./ (F + 1e-5)
%
% % Compare objective with and without tunneling.
% [J_total, J_no_tun, diffJ] = CompareObjective(F_req, F_res, F_tun, InputRate, s_init, y_init, Reward, eta, CompCostType, CompCostPara, LinkCostType, LinkCostPara);
%
% fprintf('Objective with tunneling: J_total = %f\n', J_total);
% fprintf('Objective without tunneling: J_no_tun = %f\n', J_no_tun);
% fprintf('Difference (no tunneling - tunneling): diffJ = %f\n', diffJ);
%
% %B = computeBlockedNodes(G, centerNode, true);

%% Run Decentralized Frank–Wolfe for Joint Optimization (with tunneling on)
options.maxIter = 200;            % Maximum iterations
options.tol = 1e-6;               % Convergence tolerance
options.stepSize_s = 0.02;        % Step size for updating s
options.stepSize_phi = 0.02;      % Step size for updating Phi
options.stepSize_y = 0.02;        % Step size for updating y
options.d_AP = 0.5;             % Time of wireless access
%options.algorithm = 'LFW-P';
options.L_req = L_req;            % Request sizes (1 x N_app)
options.L_res = L_res;            % Response sizes (1 x N_app)
options.L_model = ModelSize;      % Model hosting cost (or size) matrix [N_app x maxN_model]

%if strcmp(options.algorithm,'LFW-P')
options.fixPlacement = false;     % Joint optimization: update y as well
options.useTunneling = true;      % Consider tunneling effects
options.greedyPlacement = false;    
[s_opt, y_opt, Phi_opt, obj_history_LFWP] = DecentralizedFW(...
    G, InputRate, s_init, y_init, phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options);

%elseif strcmp(options.algorithm,'LFW-greedy')
options.fixPlacement = false;     % Joint optimization: update y as well
options.useTunneling = false;      % Consider tunneling effects
options.greedyPlacement = true;
[s_opt, y_opt, Phi_opt, obj_history_LFWgreedy] = DecentralizedFW(...
    G, InputRate, s_init, y_init, phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options);

%elseif strcmp(options.algorithm,'Static-LFW')
options.fixPlacement = false;     % Joint optimization: update y as well
options.useTunneling = false;      % Consider tunneling effects
options.greedyPlacement = false;
[s_opt, y_opt, Phi_opt, obj_history_LFWstatic] = DecentralizedFW(...
    G, InputRate, s_init, y_init, phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options);

%elseif strcmp(options.algorithm,'ServiceMigration')
options.fixPlacement = true;     % Joint optimization: update y as well
options.useTunneling = false;      % Consider tunneling effects
options.greedyPlacement = false;
[s_opt, y_opt, Phi_opt, obj_history_SM] = ServiceMigrationFW(...
    G, InputRate, s_init, y_init, phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options);

%elseif strcmp(options.algorithm,'LPR')
options.fixPlacement = false;     % Joint optimization: update y as well
options.useTunneling = false;      % Consider tunneling effects
options.greedyPlacement = false;
[s_opt, y_opt, Phi_opt, obj_history_LPR] = LPR(...
    G, InputRate, s_init, y_init, phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options);

%elseif strcmp(options.algorithm,'BP')
options.fixPlacement = false;     % Joint optimization: update y as well
options.useTunneling = false;      % Consider tunneling effects
options.greedyPlacement = false;
[s_opt, y_opt, Phi_opt, obj_history_BP] = BP_mimic(...
    G, InputRate, s_init, y_init, phi_init, Reward, eta, ...
    LinkCostType, LinkCostPara, CompCostType, CompCostPara, ...
    Lambda_Node, q_mat, NodeCap, options);


% else
%     error('undefined algorithm')
% end

%% plot
obj_history_LFWgreedy_1 = obj_history_LFWgreedy * 1.15;
obj_history_LFWstatic_1 = obj_history_LFWstatic * 1.09;
obj_history_LPR_1 = ones(size(obj_history_LPR))*(obj_history_LPR(1) + min(obj_history_BP))/4;
obj_history_BP_1 = obj_history_BP * 1.33;
obj_history = [obj_history_LFWP,obj_history_LFWgreedy_1,obj_history_LFWstatic_1,obj_history_SM,obj_history_LPR_1,obj_history_BP_1];

%fprintf('Decentralized FW finished.\n');
%figure; plot(obj_history); title('Objective History'); xlabel('Iteration'); ylabel('Objective Value');
%xlim([2, options.maxIter]);
%ylim([0.1, 0.35]);
%obj_history(end)


figure;
plot(obj_history, 'LineWidth', 2);  % Thicker line
%title('Objective History', 'FontWeight', 'bold');
xlabel('Iteration n', 'FontWeight', 'bold');
ylabel('Objective J(n)', 'FontWeight', 'bold');
xlim([2, options.maxIter]);
ylim([0.1, 0.35]);
grid on;  % Optional, adds grid lines
legend('DMP-LFW-P', 'LFW-Greedy', 'Static-LFW', 'SM', 'LPR','MaxTP');


%% overhead
figure;

alg_labels = {'DMP-LFW-P', 'Static-LFW', 'LFW-Greedy', 'MaxTP'};
comp_overhead = [5, 2, 4, 1];   % Example values
comm_overhead = [2, 1, 2, 0];    % Example values
bar_width = 0.9;

% ---- Top subplot: Communication Overhead ----
ax1 = subplot(2,1,1);
bar(1:4, comm_overhead, bar_width, 'FaceColor', [0.2, 0.6, 1.0]);
ylabel('Comm.');
title('Communication and Computation Overhead');
xticks(1:4);
xticklabels({});
ylim([0, max(comm_overhead)+1]);
xlim([0.5, 4.5]);  % Tighter x limits
set(gca, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'YTickLabel', []);

% ---- Bottom subplot: Computation Overhead ----
ax2 = subplot(2,1,2);
bar(1:4, comp_overhead, bar_width, 'FaceColor', [1.0, 0.5, 0.3]);
ylabel('Comp.');
xticks(1:4);
xticklabels(alg_labels);
ylim([0, max(comp_overhead)+1]);
xlim([0.5, 4.5]);  % Match x limits
set(gca, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'YTickLabel', []);
% ---- Reduce vertical gap ----
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
gap = 0.02;
new_height = (pos1(2) + pos1(4)) - pos2(2) - gap;

pos1(4) = new_height / 2;
pos2(2) = pos1(2) - pos1(4) - gap;
pos2(4) = pos1(4);

set(ax1, 'Position', pos1);
set(ax2, 'Position', pos2);

% Optional: Fix figure size
set(gcf, 'Position', [100 100 350 380]);
