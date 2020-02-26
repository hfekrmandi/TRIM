function simple_distrib_information_atmosphere()

clear all
global opt_dist
close all
dbstop if error
problem_def();
flag_converged = 0;
for i_step = 1:60
    i_step
    opt_dist.i_step = i_step;
%     flag_converged = 0;
    sim_system();
    pred();
    consenus();
    update_();
    post_process()
end
save('example_agents_CI_Comp')
%     while ~flag_converged;
%         
%         
%         % information update
%         i_time = (opt_dist.i_step-1)*(opt_dist.nIterations)+1;
%         opt_dist.i_time = i_time;
%         dataG = make_prediction(dataG);
%         
%         for i_time = (opt_dist.i_step-1)*(opt_dist.nIterations)+1:(opt_dist.i_step)*(opt_dist.nIterations)
%             dataG.x_gt(:,i_time)  = dataG.x_gt(:,(opt_dist.i_step-1)*(opt_dist.nIterations)+1);
%             opt_dist.i_time = i_time;
%             dataG.switch=0;
%             %             if i_time == 11 || i_time == 29
%             %                 dataG.switch=1;
%             %             end
%             %             if i_time > 10 & i_time<30
%             %                 dataG.Graph = opt_dist.Graphs.G_fault;
%             %             end
%             %             if i_time>30
%             %                 dataG.Graph = opt_dist.Graphs.G;
%             %             end
%             %             i_time
%             %             opt_dist.x_gt(:,i_time) =f(opt_dist.x_gt(:,max(i_time-1,1)),opt_dist.P_est(:,:,1,max(i_time-1,1)),1);
%             dataG = calc_update(dataG);
%             dataG = update_state_iterative(dataG);
%             %             plot_error_diagrams(dataG);
%         end
%         flag_converged = 1;
%     end
%    
%     for i_agent=1:opt_dist.nAgents
%         [dataG.x_est(:,i_agent,max((opt_dist.i_step)*(opt_dist.nIterations),1)),dataG.P_est(:,:,i_agent,max((opt_dist.i_step)*(opt_dist.nIterations)),1)] =  f(dataG.x_est(:,i_agent,i_time),dataG.P_est(:,:,i_agent,i_time),0);
%     end
%         [dataG.x_est_cen(:,(opt_dist.i_step)*(opt_dist.nIterations)),dataG.P_est_cen(:,:,(opt_dist.i_step)*(opt_dist.nIterations))]=...
%                   f(dataG.x_est_cen(:,(opt_dist.i_step)*(opt_dist.nIterations)),dataG.P_est_cen(:,:,(opt_dist.i_step)*(opt_dist.nIterations)),1);
%     [dataG.x_gt(:,(opt_dist.i_step)*(opt_dist.nIterations)+1),~]=f(dataG.x_gt(:,(opt_dist.i_step)*(opt_dist.nIterations)),dataG.P_est_cen(:,:,(opt_dist.i_step)*(opt_dist.nIterations)),1);
%      plot_error_diagrams(dataG);
%     plot_result(dataG.P_est,dataG.x_est,dataG.e_est)
%     
%     
% end

end



















