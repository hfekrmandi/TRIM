function consenus()
global opt_dist
% initialize the concsensus values
[opt_dist.result.pred.x_cen,opt_dist.result.pred.P_cen] = f(opt_dist.result.prior.x_cen,opt_dist.result.prior.P_cen,1);
opt_dist.result.pred.x_cen = opt_dist.result.prior.x_cen;
opt_dist.result.pred.P_cen = opt_dist.result.prior.P_cen;

opt_dist.result.pred.Y_cen = inv(opt_dist.result.pred.P_cen);
opt_dist.result.pred.y_cen = opt_dist.result.pred.Y_cen*opt_dist.result.pred.x_cen;

for i_agent = 1 : opt_dist.nAgents
    [opt_dist.result.pred.x_bar(:,i_agent),opt_dist.result.pred.P_bar(:,:,i_agent)] = ...
        f(opt_dist.result.prior.x_bar(:,i_agent),opt_dist.result.prior.P_bar(:,:,i_agent),1);
    
    opt_dist.result.pred.Y_bar(:,:,i_agent) = pinv(opt_dist.result.pred.P_bar(:,:,i_agent));
    opt_dist.result.pred.y_bar(:,i_agent) = opt_dist.result.pred.Y_bar(:,:,i_agent)*opt_dist.result.pred.x_bar(:,i_agent);
end