function update_()
global opt_dist
% centralized update
opt_dist.result.update.delta_i_cen = zeros(size(opt_dist.result.pred.x_cen));
opt_dist.result.update.delta_I_cen  = zeros(size(opt_dist.result.pred.P_cen));

% for i_agent = 1 : opt_dist.nAgents
%     opt_dist.result.update.delta_i_cen =opt_dist.result.update.delta_i_cen  +  opt_dist.result.consenus{i_agent}.delta_i{1};
%     opt_dist.result.update.delta_I_cen =opt_dist.result.update.delta_I_cen  +  opt_dist.result.consenus{i_agent}.delta_I{1};
% end
% H_cen = opt_dist.result.obs.H_cen;
% z_cen = opt_dist.sim.obs.z_cen;
% r_var_cen =    opt_dist.sim.obs.r_var_cen;
[av_delta_I,av_delta_i] = calculate_aver_info(opt_dist.result);


opt_dist.result.update.delta_I_cen = 9*av_delta_I;
opt_dist.result.update.delta_i_cen = 9*av_delta_i;



% decentralized update
for i_consensus=1:opt_dist.nSteps ;
    opt_dist.result.est{i_consensus}.Y_cen = opt_dist.result.pred.Y_cen + opt_dist.result.update.delta_I_cen;
    opt_dist.result.est{i_consensus}.y_cen = opt_dist.result.pred.y_cen + opt_dist.result.update.delta_i_cen;
    
    for j_agent = 1 : opt_dist.nAgents
        if opt_dist.FLAGS.our_method
            opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent) = opt_dist.result.consenus{j_agent}.Y_prior{i_consensus}+...
                opt_dist.nAgents*opt_dist.result.consenus{j_agent}.delta_I{i_consensus};
            opt_dist.result.est{i_consensus}.y_bar(:,j_agent) = opt_dist.result.consenus{j_agent}.y_prior{i_consensus}+...
                opt_dist.nAgents*opt_dist.result.consenus{j_agent}.delta_i{i_consensus};
        end
        if opt_dist.FLAGS.compare_with_CI
            opt_dist.result.est{i_consensus}.Y_bar_CI(:,:,j_agent) = opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus};
            opt_dist.result.est{i_consensus}.y_bar_CI(:,j_agent) = opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus};
        end
    end
end

if opt_dist.FLAGS.debug
    disp(['RCOND P_cen prediction = ',num2str(rcond(opt_dist.result.est{i_consensus}.Y_cen ))])
    if rcond(opt_dist.result.est{i_consensus}.Y_cen )<=0.0001
        figure(opt_dist.figures.fig_cov_debug)
    end
end



% opt_dist.result.est.x_cen

% decen update
end