function dataG = make_prediction(dataG)
global opt_dist
i_time = opt_dist.i_time;
%%
% [dataG.x_pred_cen(:,i_time),dataG.P_pred_cen(:,:,i_time)] = f(dataG.x_est_cen(:,max(i_time-1,1)),dataG.P_est_cen(:,:,max(i_time-1,1)),1);
% dataG.Y_pred_cen(:,:,i_time) = pinv(dataG.P_pred_cen(:,:,i_time));
% dataG.y_pred_cen(:,i_time) = dataG.Y_pred_cen(:,:,i_time)*dataG.x_pred_cen(:,i_time);

dataG.x_pred_cen(:,i_time) = dataG.x_est_cen(:,max(i_time-1,1));
dataG.P_pred_cen(:,:,i_time) = dataG.P_est_cen(:,:,max(i_time-1,1));

dataG.Y_pred_cen(:,:,i_time) = pinv(dataG.P_pred_cen(:,:,i_time));
dataG.y_pred_cen(:,i_time) = dataG.Y_pred_cen(:,:,i_time)*dataG.x_pred_cen(:,i_time);




for i_agent = 1 : opt_dist.nAgents
    %     [dataG.x_pred(:,i_agent,i_time),dataG.P_pred(:,:,i_agent,i_time)] = f(dataG.x_est(:,i_agent,max(i_time-1,1)),dataG.P_est(:,:,i_agent,max(i_time-1,1)),1);
    dataG.x_pred(:,i_agent,i_time) = dataG.x_est(:,i_agent,max(i_time-1,1));
    dataG.P_pred(:,:,i_agent,i_time) = dataG.P_est(:,:,i_agent,max(i_time-1,1));
    
    dataG.Y_pred(:,:,i_agent,i_time) = pinv(dataG.P_pred(:,:,i_agent,i_time));
    dataG.y_pred(:,i_agent,i_time) = dataG.Y_pred(:,:,i_agent,i_time)*dataG.x_pred(:,i_agent,i_time);
end
end