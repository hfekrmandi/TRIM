function agents_state_plot(location)
% Inputs: 
    % location - String for the matrix of workspace data to import
% Outputs: 
    % None
% Description:
    % Plots the agent states vs. time
        % Figure 1 compares the state values between different methods
        % Figure 2 compares the state errors (difference between each method 
        % and the ground truth) between different methods
    % In each figure:
        % gt is the ground truth (True values)
        % gold is the centralized estimates
        % hybrid is the hybrid method estimates
        % ici is the ICI method estimates

%%
% clearvars
% close all
% clc

% Load in the data from a prior run
% for cubesat testing, this is 'cubesat_test_3.mat'
if nargin == 0
    location = 'cubesat_test_3.mat';
end
load(location)

% Create the number of agents and the number of time steps data
nAgents = opt_dist.nAgents;
nSteps = opt_dist.nSteps;
fin_step = opt_dist.i_step;
agents = 1:nAgents;
time = 1:fin_step;

% Initialize the data arrays
T = size(x_gt);
% gt is the ground truth (True values)
gt_data = zeros(fin_step, nAgents);
% gold is the centralized estimates
gold_data = zeros(fin_step, nAgents);
% hybrid is the hybrid method estimates
hyb_data = zeros(fin_step, nAgents);
% ici is the ICI method estimates
ici_data = zeros(fin_step, nAgents);

% Fill the data arrays with stored data
for i=1:fin_step
    gt_data(i,:) = x_gt{:,i};
    gold_data(i,:) = x_gold{i,:};
    hyb_data(i,:) = x_Hybrid{i,:};
    ici_data(i,:) = x_ICI{i,:};
end

%% Create figure 1 (State values)
fig = figure(1);
hold on
colors = ['k', 'b', 'r', 'g'];
for i = 1:nAgents
    plot(time, gt_data(:,i)', 'DisplayName', strcat('Ground Truth Agent ', int2str(i)));%, colors(1));
    plot(time, gold_data(:,i), 'DisplayName', strcat('Centralized Agent ', int2str(i)));%, colors(2));
    plot(time, hyb_data(:,i), 'DisplayName', strcat('Distributed Agent ', int2str(i)));%, colors(3));
    %plot(time, ici_data(:,i), colors(4));
end
legend();
xlabel('Time Steps');
ylabel('State Values');

%% Create figure 2 (State errors)
fig = figure(2);
hold on
temp = (gt_data(:,:)) - gold_data(:,:);
for i = 1:nAgents
    plot(time, gt_data(:,i) - gold_data(:,i), 'DisplayName', strcat('Centralized Error Agent ', int2str(i)));%, colors(2));
    plot(time, gt_data(:,i) - hyb_data(:,i), 'DisplayName', strcat('Distributed Error Agent ', int2str(i)));%, colors(3));
    %plot(time, gt_data(:,i) - ici_data(:,i), colors(4));
end

y_max = max([gt_data - gold_data, gt_data - hyb_data], [], 'all');
y_min = min([gt_data - gold_data, gt_data - hyb_data], [], 'all');
x_vals = [[2, 2, 4, 4]; [6, 6, 8, 8]; [100 100 100 100]] - 0.5;
j_max = 2;
j = 1;
for i = 1:nSteps
    if j <= j_max && (x_vals(j,1) <= i && x_vals(j+1,1) >= i)
        x = x_vals(j,:);
        j = j + 1;
        y = [y_min, y_max, y_max, y_min];
        fill(x, y, 'r', 'facealpha', 0.1,'edgecolor','none','HandleVisibility','off');
    end
end
legend('Location','northeastoutside');
xlabel('Time Steps');
ylabel('State Value Error (Absolute)');

%% Create Figure 3 (Consensus errors)
figure(3);
P_hyb = zeros([nAgents, nAgents, nSteps]);
x_hyb_err = zeros([nAgents, nSteps]);
P_ici = zeros([nAgents, nAgents, nSteps]);
x_ici_err = zeros([nAgents, nSteps]);
steps = 1:(nSteps*fin_step);
for i = 1:nAgents
    for j = 1:fin_step
        start = nSteps * (j - 1);
        for k = 1:nSteps
            Y_hyb = consensus_history{j}{i}.Y_prior{k} + ...
                       consensus_history{j}{i}.delta_I{k};
            y_hyb = consensus_history{j}{i}.y_prior{k} + ...
                       consensus_history{j}{i}.delta_i{k};
            P_hyb(:,:,k+start) = pinv(Y_hyb);
            x_hyb_err(:,k+start) = x_gt{j} - P_hyb(:,:,k+start)*y_hyb;

            Y_ici = consensus_history{j}{i}.Y_prior_CI{k} + ...
                       consensus_history{j}{i}.delta_I_CI{k};
            y_ici = consensus_history{j}{i}.y_prior_CI{k} + ...
                       consensus_history{j}{i}.delta_i_CI{k};
            P_ici(:,:,k+start) = pinv(Y_ici);
            x_ici_err(:,k+start) = x_gt{j} - P_ici(:,:,k+start)*y_ici;
        end
    end
    %subplot(nAgents, 1, i);
    ylabel('State Value Error (Absolute)');
    hold on
    plot(steps, x_hyb_err(1,:), 'DisplayName', strcat('Agent ', int2str(i), ' Estimating Agent ', '1'));%, colors(3));
    plot(steps, x_hyb_err(2,:), 'DisplayName', strcat('Agent ', int2str(i), ' Estimating Agent ', '2'));%, colors(3));
    %plot(steps, x_ici_err, colors(4));
    y_max = max([x_hyb_err, x_ici_err], [], 'all');
    y_min = min([x_hyb_err, x_ici_err], [], 'all');
    if fin_step > 1
        x = [0, 0, nSteps, nSteps];
        y = [y_min, y_max, y_max, y_min];
        for j = 1:ceil(fin_step / 2)
            fill(x, y, 'g', 'facealpha', 0.2,'edgecolor','none','HandleVisibility','off');
            x = x + 2*nSteps;
        end
    end
end
legend('Location','northeastoutside');
xlabel('Consensus Steps');
ylabel('Distributed State Value Error vs. Consensus steps');

end