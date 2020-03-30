function [perf_index_prob,perf_index] = gold_standard()
addpath(genpath(pwd));
clear all
global opt_dist
profile on
close all
% dbclear all
dbstop if error
% dbstop if warning
% Fix the random number generator seed to make the runs repeatable
% Feel free to change this.
rng(0)

% Major Steps - 
% 1) create_sys_gold
    % Creates a random, stable system
% 2) problem_def_gold
    % Defines first-time variables and the graph
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
        % result.est.Y_cen, result.est{}.y_cen
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
            % Type: Float
    % Graphs
        % Type: Struct of GraphsGraph
    % x_gt
        % Type: Matrix (nx1)
    % result
        % gt
            % x_bar
                % Type: Matrix (nx1)
        % prior
            % x_cen
                % Type: Matrix (nx1)
            % P_cen
                % Type: Matrix (nxn)
            % x_bar
                % Type: Matrix (nxp)
            % P_bar
                % Type: Matrix (nxn xp)
            % x_bar_CI
                % Type: Matrix (nxp)
            % P_bar_CI
                % Type: Matrix (nxn xp)
        % initial
            % x_bar
                % Type: Matrix (nxp)
            % P_bar
                % Type: Matrix (nxn xp)
            % x_bar_CI
                % Type: Matrix (nxp)
            % P_bar_CI
                % Type: Matrix (nxn xp)
        % pred
            % x_cen
                % Type: Matrix (nx1)
            % P_cen
                % Type: Matrix (nxn)
            % Y_cen
                % Type: Matrix (nxn)
            % y_cen
                % Type: Matrix (nx1)
        % obs
            % h_cen
                % Type: Matrix (px1)
            % H_cen
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
    
% This is the size of the reduced atmospheric system
%reduced_dimention_size = 2;
%[A,x0,B,C] = create_sys_atmosphere_gold(reduced_dimention_size);

% create_sys_gold(n,p,m) generates an n-th order model with p outputs and m inputs.
[A,x0,B,C] = create_sys_gold(3,5,1);


flag_converged = 0;
global fail_prob reg_deg
% range_prob = [ 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
% range_prob = [ 0.4 0.6 0.8  1];
% range_prob = [ 0 0.2  0.4 0.6 0.8 1];
range_prob = [ 0.2 ];

% range_prob = [[0:0.2:0.4],[0.5:0.05:0.8],0.9,1];

% range_reg = [ 2 4 6 8];
range_reg = [ 4];

% range_prob = [ 1];
problem_def_gold(A,B,C,x0);

if strcmp(opt_dist.scenario, '1')
    for i_step = 1:5
        i_step
        opt_dist.i_step = i_step;
        %     flag_converged = 0;
        sim_system_gold();
        pred();
        consenus_gold();
        calc_super_gold_update();
        [error_results{i_step}] = post_process_gold();
    end
    profile viewer
    
else
    range_step =24;
    for j_reg=1:length(range_reg)
        reg_deg = range_reg(j_reg);
        tic
        for i_prob=1:length(range_prob)
            if ~(j_reg==1 && i_prob==1)
                fields = {'obs' ,'result','Graph_History','sim'};
                opt_dist = rmfield(opt_dist,fields);
            end
            problem_def_gold(A,B,C,x0);
            
            fail_prob =  range_prob(i_prob);
            opt_dist.reg_degree = reg_deg;
            
            for i_step = 1:range_step
                i_step
                opt_dist.i_step = i_step;
                sim_system_gold();
                pred();
                consenus_gold();
                calc_super_gold_update();
                time_(j_reg,i_prob,i_step) = toc;
                [error_results{j_reg,i_prob,i_step}] = post_process_gold2();
                if (error_results{j_reg,i_prob,i_step}.error_Hybrid.e_BC_dist_cent - error_results{j_reg,i_prob,i_step}.error_Hybrid.e_BC_dist_gold_vs_cent)> 0.001
                    disp('check')
                end
            end
        end
        
    end
end
mean_ = calc_composite_results_gold(error_results,length(range_reg),length(range_prob),range_step)

end

