function [perf_index_prob,perf_index] = gold_standard()
addpath(genpath(pwd));
clear all
profile on
close all
dbstop if error
global opt_dist

% This is the size of the reduced atmospheric system
error_index = 1; %initialising index for error_ variable.

conv_index = 1;

% Fix the random number generator seed to make the runs repeatable
% Feel free to change this.
rng(0)

% create_sys_gold(n,p,m) generates an n-th order model with p outputs and m inputs.
[A,x0,B,C] = create_sys_gold(9,9,3);
disp('controllability:'); disp(rank(ctrb(A,B)));
disp('observability:'); disp(rank(obsv(A,C)));
%% Problem posed
global fail_prob reg_deg

range_prob = [ 0.02 ]; %0.2
range_reg = [ 4];
converg_steps = 60;
fail_prob = 0.2;
problem_def_gold(A,B,C,x0,converg_steps);

for i_step = 1:5
    i_step
    opt_dist.i_step = i_step;

    sim_system_gold();

    pred();

    consenus_gold();
    calc_super_gold_update();
    [error_results{i_step}] = post_process_gold();
end
profile viewer


[error_(error_index,conv_index),mean_(error_index,conv_index)] = calc_composite_results_gold(error_results,length(range_reg),length(range_prob),range_step);
time_array(error_index,conv_index) = mean(squeeze(time_));

assignin('base','mean_',mean_);%store variable in workspace
assignin('base','error_',error_);
assignin('base','n_receptors_array',n_receptors_array);
assignin('base','converg_steps_array',converg_steps_array);
assignin('base','time_array',time_array);
assignin('base','P_gold',P_gold);
assignin('base','P_ICI',P_ICI);
assignin('base','P_Hybrid',P_Hybrid);
assignin('base','x_gold',x_gold);
assignin('base','x_ICI',x_ICI);
assignin('base','x_Hybrid',x_Hybrid);

end
