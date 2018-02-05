function atmosphere_for_paper()
clear all
global opt_dist
close all
dbstop if error
dbstop if warning

problem_def();

flag_converged = 0;
for i_step = 1:70
    i_step
    opt_dist.i_step = i_step;
%     flag_converged = 0;
    sim_system();
    pred();
    consenus();
    update_();
    post_process()
end
save('example_agents_CI_Comp_70_steps_4')
end
