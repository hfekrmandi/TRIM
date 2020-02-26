function dataG = update_state_iterative(dataG)
global opt_dist
i_time = opt_dist.i_time;


dataG.mean_i_update = mean( dataG.i_update(:,:,i_time)' )';
dataG.mean_I_update = mean( dataG.I_update(:,:,:,i_time),3  );

dataG.Y_est_cen(:,:,i_time)=  dataG.Y_pred_cen(:,:,(opt_dist.i_step-1)*(opt_dist.nIterations)+1) + opt_dist.nAgents*dataG.mean_I_update;
dataG.y_est_cen(:,i_time)=  dataG.y_pred_cen(:,(opt_dist.i_step-1)*(opt_dist.nIterations)+1) + opt_dist.nAgents*dataG.mean_i_update;
dataG.P_est_cen(:,:,i_time)= pinv( dataG.Y_est_cen(:,:,i_time));
dataG.x_est_cen(:,i_time)=  dataG.P_est_cen(:,:,i_time)* dataG.y_est_cen(:,i_time);
dataG.e_est_cen(:,i_time) =  dataG.x_est_cen(:,i_time) -  dataG.x_gt(:,(opt_dist.i_step-1)*(opt_dist.nIterations)+1);

dataG.e_est_norm_cen(i_time) = norm( dataG.e_est_cen(:,i_time));

for i_agent = 1 : opt_dist.nAgents
    dataG.I_con(:,:,i_agent,i_time) = zeros(size(squeeze(dataG.I_update(:,:,i_agent,i_time))));
    dataG.i_con(:,i_agent,i_time) = zeros(size(squeeze(dataG.i_update(:,i_agent,i_time))));
    for  j_agent = 1 : opt_dist.nAgents
        dataG.I_con(:,:,i_agent,i_time) = dataG.I_con(:,:,i_agent,i_time) + dataG.Graph.p(i_agent,j_agent)*dataG.I_con(:,:,j_agent,max(i_time-1,1));
        dataG.i_con(:,i_agent,i_time) = dataG.i_con(:,i_agent,i_time) + dataG.Graph.p(i_agent,j_agent)*dataG.i_con(:,j_agent,max(i_time-1,1));
    end
    %%
    dataG.Y_est(:,:,i_agent,i_time)=  dataG.Y_pred(:,:,i_agent,(opt_dist.i_step-1)*(opt_dist.nIterations)+1) + opt_dist.nAgents*dataG.I_con(:,:,i_agent,i_time);
    dataG.y_est(:,i_agent,i_time)=  dataG.y_pred(:,i_agent,(opt_dist.i_step-1)*(opt_dist.nIterations)+1) + opt_dist.nAgents*dataG.i_con(:,i_agent,i_time);
    dataG.P_est(:,:,i_agent,i_time)= pinv( dataG.Y_est(:,:,i_agent,i_time));
    dataG.x_est(:,i_agent,i_time)=  dataG.P_est(:,:,i_agent,i_time)* dataG.y_est(:,i_agent,i_time);
    dataG.e_est(:,i_agent,i_time) =  dataG.x_est(:,i_agent,i_time) -  dataG.x_gt(:,(opt_dist.i_step-1)*(opt_dist.nIterations)+1);
    dataG.e_est_norm(i_agent,i_time) = norm( dataG.e_est(:,i_agent,i_time));
    %%000000000000000000000000000000000000000000000000000000000000000000000000000000
    dataG.e_i_con(i_agent,i_time) = norm( dataG.i_con(:,i_agent,i_time) -  dataG.mean_i_update);
    dataG.e_I_con(i_agent,i_time) = norm( dataG.I_con(:,:,i_agent,i_time) -  dataG.mean_I_update);
    dataG.e_P_est(:,:,i_agent,i_time) = dataG.P_est(:,:,i_agent,i_time) - dataG.P_est_cen(:,:,i_time);
    dataG.e_P_est_norm(i_agent,i_time) = norm(dataG.e_P_est(:,:,i_agent,i_time));

end
%%
end