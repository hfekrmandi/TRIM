function dataG = make_prediction_initial()
global opt_dist
% opt_dist.Graphs.G,opt_dist.x_gt
dataG.Graph = opt_dist.Graphs.G;

i_time = opt_dist.i_time;
dataG.x_gt = opt_dist.x_gt;







dataG.x_temp_cen  = dataG.x_gt;
dataG.p_temp_cen =  diag(0.05*ones(opt_dist.dimState,1));
%         [x_pred(:,i_agent,i_time),P_pred(:,:,i_agent,i_time)] = f(x_est(:,i_agent,1),P_est(:,:,i_agent,1),1);
dataG.x_pred_cen(:,i_time) = dataG.x_temp_cen;
dataG.P_pred_cen(:,:,i_time) = dataG.p_temp_cen;

dataG.Y_pred_cen(:,:,i_time) = pinv(dataG.P_pred_cen(:,:,i_time));
dataG.y_pred_cen(:,i_time) = dataG.Y_pred_cen(:,:,i_time)*dataG.x_pred_cen(:,i_time);

dataG.x_est_cen(:,i_time) = dataG.x_temp_cen;
dataG.P_est_cen(:,:,i_time) = dataG.p_temp_cen;

dataG.Y_est_cen(:,:,i_time) = pinv(dataG.P_est_cen(:,:,i_time));
dataG.y_est_cen(:,i_time) = dataG.Y_est_cen(:,:,i_time)*dataG.x_est_cen(:,i_time);

bigUncertainty = 3;
for i_agent = 1 : opt_dist.nAgents
    idx_ith_agent = zeros(opt_dist.dimState,1);
    idx_ith_agent((i_agent-1)*opt_dist.dimAgents + 1: (i_agent)*opt_dist.dimAgents) = 1;
%     dataG.x_temp  = dataG.x_gt.*idx_ith_agent;
%     dataG.x_temp  = dataG.x_gt + randn(size(dataG.x_gt,1),1)*opt_dist.motion.Q;

    P_init = bigUncertainty*diag(~idx_ith_agent) + diag(0.05.*idx_ith_agent);
    dataG.x_temp = [];
    dataG.p_temp = [];
            [dataG.x_temp ,dataG.p_temp] = f(dataG.x_gt,dataG.p_temp_cen,1);
    dataG.x_pred(:,i_agent,i_time) = dataG.x_temp_cen;
    dataG.P_pred(:,:,i_agent,i_time) = dataG.p_temp_cen;
    dataG.x_est(:,i_agent,i_time) = dataG.x_temp;
    dataG.P_est(:,:,i_agent,i_time) = dataG.p_temp;
    dataG.Y_pred(:,:,i_agent,i_time) = pinv(dataG.P_pred(:,:,i_agent,i_time));
    dataG.y_pred(:,i_agent,i_time) = dataG.Y_pred(:,:,i_agent,i_time)*dataG.x_pred(:,i_agent,i_time);
end
% dataG.h_plot = figure;
% dataG.h_error = figure;
% dataG.h_error2 = figure;

end
