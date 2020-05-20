function [perf_index_prob,perf_index] = gold_standard()
clear all
global opt_dist
profile on
close all
% dbclear all
dbstop if error
% dbstop if warning
%rng(0)


% Major Steps - 
% 1) create_sys_gold
    % Creates a random, stable system
% 2) problem_def_gold
    % Defines first-time variables and the connection graph
% Loop through a number of time steps:
% 3) sim_system_gold();
    % Simulates a time-step of the system
    % a) f_sim
        % computes x[k+1] = A*x[k] + B*u[k] (u is Q)
    % b) sim_obs_gold
        % assigns z, r_var, and id_vis
% 4) pred();
    % Computes x, P, Y, h, H, and id_vis for both centralized (_cen) and HCCI (_bar)
    % Paper steps: EQ 2/3, EQ 8b
    % a) f
        % computes x[k+1] = A*x[k] + B*u[k] (u is Q) and
        % P[k+1] = A*P[k]*A' + u[k]*I
    % b) h_calc_ver2
        % computes centralized h, H, id_vis
    % c) h_calc_gold
        % computes HCCI h, H, id_vis
% 5) consenus_gold();
    % The main consensus step
    % a) calc_ci_weights_ver3
        % 
% 6) calc_super_gold_update();
    % Plots data
    % a) f
        % computes x[k+1] = A*x[k] + B*u[k] (u is Q) and
        % P[k+1] = A*P[k]*A' + u[k]*I

% function input/output and in-depth descriptions:
% function [A,x0,B,C] = create_sys_gold(n,p,m)
    % Description:
        % Creates a random n-th order system with p outputs and m inputs.
        % Makes sure the system is stable, but there is no guaranty on
        % controllability/observability
    % Inputs:
        % n - The order of the system
        % p - The number of outputs of the system
            % MUST BE >= 5
        % m - The number of inputs of the system
    % Outputs:
        % A - The system matrix A (nxn)
        % x0 - The initial system state (all zeros, nx1)
        % B - The input matrix B (nxm)
        % C - The output matrix C (pxn)
        % (opt) Q - The initial input (all ones, mx1)

% function problem_def_gold(A,B,C,x0)
    % Description:
        % Defines first-time variables and the graph
        % Opens log file
    % Inputs:
        % A - The system matrix A (nxn)
        % x0 - The initial system state (all zeros, nx1)
        % B - The input matrix B (nxm)
        % C - The output matrix C (pxn)
    % Outputs:
        % (opt) Sets up each agent's initial and prior x_bar and P_bar

% function sim_system_gold()
    % Description:
        % Simulates a time-step of the system
    % All inputs/outputs are from opt_dist
    % Inputs:
        % A, B, Q, x_current
    % Outputs:
        % x_next
        % z, r_var, id_vis

% function pred()
    % Description:
        % Computes x, P, Y, h, H, and id_vis for both centralized (_cen) and HCCI (_bar)
    % All inputs/outputs are from opt_dist
    % Inputs:
        % result.prior.x_cen, result.prior.P_cen
        % sys.A, sys.B, source.Q, result.prior.x_bar, result.prior.P_bar
        % result.pred.x_bar, result.pred.P_bar, result.pred.Y_bar
        % result.pred.x_bar_CI, result.pred.P_bar_CI, result.pred.Y_bar_CI
        % result.pred.x_cen
    % Outputs:
        % result.pred.Y_cen, result.pred.x_cen, result.pred.P_cen
        % result.pred.y_cen
        % result.obs.h_cen, result.obs.H_cen, result.obs.id_vis_cen

% function consenus_gold()
    % Description:
        % The main consensus step
    % All inputs/outputs are from opt_dist
    % Inputs:
        % sim.obs.r_var, result.obs.H, result.obs.z
        % result.pred.Y_bar, result.pred.y_bar
        % opt_dist.result.consenus.Y_prior,opt_dist.result.consenus.y_prior
        % result.consenus.delta_I, result.consenus.delta_i
    % Outputs:
        % result.est.Y_cen, result.est.y_cen
        % opt_dist.result.consenus.Y_prior,opt_dist.result.consenus.y_prior
        % result.consenus.delta_I, result.consenus.delta_i
        % result.est.Y_bar, result.est.y_bar

% function calc_super_gold_update()
    % Description:
        % Stores and plots data
    % All inputs/outputs are from opt_dist
    % Inputs:
        % result.obs.H,result.obs.z, result.obs.r_var
    % Outputs:
        % opt_dist.result.est_gold.Y_bar
        % opt_dist.result.est_gold.y_bar
        % opt_dist.result.est_gold.P_bar
        % opt_dist.result.est_gold.x_bar

% function post_process_gold()
    % Description:
    % All inputs/outputs are from opt_dist
    % Inputs:
        % result.est.Y_bar, result.est.y_bar
        % result.est_gold.P_bar
        % result.est_gold.x_bar
    % Outputs:

% opt_dist:
    % sys
        % 5x1 State Space
        % The state matrices A, B, C, and D
    % A
        % Type: Matrix (nxn)
        % Value: Equivalent to sys.A
    % B
        % Type: Matrix (nxm)
        % Value: Equivalent to sys.B
    % C
        % Type: Matrix (pxn)
        % Value: Equivalent to sys.C
    % source
        % Q
            % Description: Input (u)
            % Type: Float
    % fid_log
        % Type: Integer
    % FLAGS
        % compare_with_CI
            % Type: Boolean
        % our_method
            % Type: Boolean
        % pure_CI
            % Type: Boolean
        % debug_CI
            % Type: Boolean
        % verbose
            % Type: Boolean
        % debug
            % Type: Boolean
        % debug_consensus
            % Type: Boolean
        % obs_noise_type
            % Type: String
            % Values: 'absolute', 'relative'
    % nAgents
        % Type: Integer
        % Value: p
    % dt
        % Type: Float
    % dimAgents
        % Type: Integer
    % obs
        % Range
            % Type: Float
        % rel_perc
            % Type: Float
            % relative noise percentage
        % R
            % Type: Float
            % absolute noise value
    % dimState
        % Type: Integer
    % reg_degree
        % Type: Integer
    % n_degree_graph
        % Type: Integer
    % figures
        % Type: Struct of Figures
    % iter_interest
        % Type: Integer
    % dimObs
        % Type: Integer
    % nIterations
        % Type: Integer
    % nSteps
        % Type: Integer
    % scenario
            % Type: String
            % Values: '1', '2'
    % motion
        % Q
            % Description: Input (u)
            % Type: Float
    % Graphs
        % Type: Struct of GraphsGraph
    % x_gt
        % Type: Matrix (nx1)
    % result
        % gt - (ground truth)
            % x_bar
                % Description: Current system true state
                % Type: Matrix (nx1)
        % prior
            % x_cen
                % Description: Most recent Centralized estimation of x
                % Type: Matrix (nx1)
            % P_cen
                % Description: Most recent Centralized estimation of P
                % Type: Matrix (nxn)
            % x_bar
                % Description: Most recent Hybrid estimation of x
                % Type: Matrix (nxp)
            % P_bar
                % Description: Most recent Hybrid estimation of P
                % Type: Matrix (nxn xp)
            % x_bar_CI
                % Description: Most recent CI estimation of x
                % Type: Matrix (nxp)
            % P_bar_CI
                % Description: Most recent CI estimation of P
                % Type: Matrix (nxn xp)
        % initial
            % x_bar
                % Description: First Hybrid estimation of x
                % Type: Matrix (nxp)
            % P_bar
                % Description: First Hybrid estimation of P
                % Type: Matrix (nxn xp)
            % x_bar_CI
                % Description: First CI estimation of x
                % Type: Matrix (nxp)
            % P_bar_CI
                % Description: First CI estimation of P
                % Type: Matrix (nxn xp)
        % pred
            % x_cen
                % Description: Centralized estimation of x
                % Type: Matrix (nx1)
            % P_cen
                % Description: Centralized estimation of P
                % Type: Matrix (nxn)
            % Y_cen
                % Type: Matrix (nxn)
            % y_cen
                % Type: Matrix (nx1)
            % x_bar
                % Description: Hybrid estimation of x
                % Type: Matrix (nxp)
            % P_bar
                % Description: Hybrid estimation of P
                % Type: Matrix (nxn xp)
            % Y_bar
                % Type: Matrix (nxn)
            % y_bar
                % Type: Matrix (nx1)
            % x_bar_CI
                % Description: CI estimation of x
                % Type: Matrix (nxp)
            % P_bar_CI
                % Description: CI estimation of P
                % Type: Matrix (nxn xp)
            % Y_bar_CI
                % Type: Matrix (nxn)
            % y_bar_CI
                % Type: Matrix (nx1)
        % obs
            % h_cen
                % Description: State estimates (same as x_bar?)
                % Type: Matrix (px1)
            % H_cen
                % Description: Matrix of measurable output variables (H)
                % Type: Matrix (pxn)
            % id_vis_cen
                % Type: Matrix (px1)
            % h
                % Type: Matrix (1xp)
            % H
                % Type: Matrix (1xp)
            % id_vis
                % Type: Matrix (1xp)
            % h_CI
                % Type: Matrix (1xp)
            % H_CI
                % Type: Matrix (1xp)
            % id_vis_CI
                % Type: Matrix (1xp)
        % update
            % delta_i_cen
                % Type: Matrix (nx1)
            % delta_I_cen
                % Type: Matrix (nx1)
        % est
            % 1x60 Array of:
                % Y_cen
                    % Type: Matrix (nx1)
                % y_cen
                    % Type: Matrix (nx1)
                % Y_bar
                    % Type: Matrix (nx1)
                % y_bar
                    % Type: Matrix (nx1)
        % consenus
            % 1xp Array of:
                % Y_prior
                    % 1x60 Array of:
                        % Type: Matrix (nxn)
                % y_prior
                    % 1x60 Array of:
                        % Type: Matrix (nx1)
                % delta_I
                    % 1x60 Array of:
                        % Type: Matrix (nxn)
                % delta_i
                    % 1x60 Array of:
                        % Type: Matrix (nx1)
                % group_set
                    % 1x60 Array of:
                        % Type: Matrix (px1)
    % i_time
        % Type: Float
    % i_step
        % Type: Integer
    % dataG
        % Type: Graph
    % Graph_History
        % Type: Matrix (1xi_step) of (pxp) Matrices
    % sim
        % Stores variables at each time step
        % gt
            % x_bar
                % Type: Matrix (nx1)
            % x_bar_history
                % Type: Matrix (nx1 xi_step)
        % obs
            % z_cen
                % Type: Matrix (px1 xi_step)
            % r_var_cen
                % Type: Matrix (px1 xi_step)
            % id_vis_cen
                % Type: Matrix (px1 xi_step)
            % z
                % Type: Matrix (i_stepx 1xp)
            % r_var
                % Type: Matrix (i_stepx 1xp)
            % id_vis
                % Type: Matrix (i_stepx 1xp)
            % z_CI
                % Type: Matrix (i_stepx 1xp)
            % r_CI_var
                % Type: Matrix (i_stepx 1xp)
            % id_vis_CI
                % Type: Matrix (i_stepx 1xp)
    % connection_history
        % Matrix (px(i_step - 1)) of Connection histories
        % Each history is a Matrix (1xp)

% Steps to getting cubesat sim:
    % Create state matrices (Block matrix AxI)
    % Figure out ranges and whether connection breaks with range disparity
    % Figure out plotting
addpath(genpath(pwd));
flag_converged = 0;
global fail_prob reg_deg
% Set probability of link failure at each step
range_prob = [0.1];
fail_prob = range_prob;

% Number of agents connected to each agent
% Degree of graph nodes
range_reg = [2];

% Creates a block diagonal matrix of n-th order models each with p outputs and m inputs
%[A,x0,B,C] = cubesat_create_sys_test(n, p, m, nAgents);
[A,x0,B,C] = cubesat_create_sys_test(1, 1, 1, 2);

% range_prob = [1];
problem_def_gold(A,x0(1:size(A,1)));
opt_dist.scenario = '1';

for i_step = 1:10
    i_step
    opt_dist.i_step = i_step;
    %     flag_converged = 0;
    sim_system_gold();
    pred();
    consenus_gold();
    % Store consensus data
    consensus_history{i_step} = opt_dist.result.consenus;
    % update
    calc_gold_update();

    [error_results{i_step}] = post_process_gold();
    for i_agent=1:opt_dist.nAgents

        x_gt{i_step} = opt_dist.result.gt.x_bar;
        
        P_gold{i_step,i_agent} = opt_dist.result.est_gold{i_agent}.P_bar;
        x_gold{i_step,i_agent} = opt_dist.result.est_gold{i_agent}.x_bar;

        P_Hybrid{i_step,i_agent} = inv(opt_dist.result.est{opt_dist.nSteps}.Y_bar(:,:,i_agent));
        x_Hybrid{i_step,i_agent} = P_Hybrid{i_step,i_agent}*(opt_dist.result.est{opt_dist.nSteps}.y_bar(:,i_agent));

        P_ICI{i_step,i_agent} = inv(opt_dist.result.est{opt_dist.nSteps}.Y_bar_CI(:,:,i_agent));
        x_ICI{i_step,i_agent} = P_ICI{i_step,i_agent}*(opt_dist.result.est{opt_dist.nSteps}.y_bar_CI(:,i_agent));

    end
end
% plot_error(error_);

assignin('base','error_results',error_results);
assignin('base','opt_dist', opt_dist);
assignin('base', 'consensus_history', consensus_history);

assignin('base','P_gold',P_gold);
assignin('base','P_ICI',P_ICI);
assignin('base','P_Hybrid',P_Hybrid);
assignin('base','x_gt',x_gt);
assignin('base','x_gold',x_gold);
assignin('base','x_ICI',x_ICI);
assignin('base','x_Hybrid',x_Hybrid);
end