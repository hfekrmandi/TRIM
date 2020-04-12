function [perf_index_prob,perf_index] = gold_standard()
addpath(genpath(pwd));
clear all
close all

profile on
dbstop if error
global opt_dist

% Fix the random number generator seed to make the runs repeatable
% Feel free to change this.
rng(0)

% create_sys_gold(n,p,m) generates an n-th order model with p outputs and m inputs.
[A,x0,B,C] = create_sys_gold(9,9,3);
disp('controllability:'); disp(rank(ctrb(A,B)));
disp('observability:'); disp(rank(obsv(A,C)));
%% Problem posed
converg_steps = 60;
problem_def_gold(A,B,C,x0,converg_steps);
range_step = 3;
global fail_prob 
range_prob = [opt_dist.reg_degree];
range_reg = [0.2];
fail_prob = range_prob(1);

for i_step = 1:range_step
    i_step
    opt_dist.i_step = i_step;
    sim_system_gold();
    pred();
    consenus_gold();
    calc_super_gold_update();
    [error_results{i_step}] = post_process_gold();
end
profile viewer
end
