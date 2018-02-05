function calc_update()
global opt_dist

for i_agent = 1 : opt_dist.nAgents
    [z,h,H,id_vis] = h_calc(dataG.x_gt(:,(opt_dist.i_step-1)*(opt_dist.nIterations)+1),dataG.x_pred(:,i_agent,(opt_dist.i_step-1)*(opt_dist.nIterations)+1),i_agent,0);
    opt_dist.result.z{i_agent} =
    opt_dist.result.h{i_agent} =
    opt_dist.result.H{i_agent} =
    opt_dist.result.id_vis{i_agent} =
    dataG.I_update(:,:,i_agent,i_time) = 1/(opt_dist.obs.R)*(H'*H);
    dataG.i_update(:,i_agent,i_time) = 1/(opt_dist.obs.R)*H'*z;
    
    dataG.I_con(:,:,i_agent,(opt_dist.i_step-1)*(opt_dist.nIterations)+1) = dataG.I_update(:,:,i_agent,i_time);
    dataG.i_con(:,i_agent,(opt_dist.i_step-1)*(opt_dist.nIterations)+1) = dataG.i_update(:,i_agent,i_time);
    
    
end
end