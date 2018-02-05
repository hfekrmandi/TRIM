function h_calc()

global opt_dist

z = [];
id_vis = [];
H = [];
h = [];

for idx_agent=1:opt_dist.nAgents
    x_pred = opt_dist.result.pred.x_bar(:,idx_agent);
    if opt_dist.FLAGS.compare_with_CI
        x_pred_CI = opt_dist.result.pred.x_bar_CI(:,idx_agent);
        id_vis_CI = [];
        H_CI = [];
        h_CI = [];
    end
    idx_obs_agent = find(opt_dist.Graphs.G_obs.Adj(idx_agent,:)~=0);
    % idx_obs{1} = [1 2];
    % idx_obs{2} = [1 2 3];
    % idx_obs{3} = [2 3 4];
    % idx_obs{4} = [3 4];
    id_vis = [];
    H = [];
    h = [];
    
    for i=idx_obs_agent
        H_temp = zeros(opt_dist.dimObs,size(opt_dist.x_gt,1));
        

            h_temp = opt_dist.C(idx_agent,:) * x_pred;
            %         if norm(z_temp) <=opt_dist.obs.Range
            id_vis = [id_vis,i];
            h = [h;h_temp];
            H_temp = opt_dist.C(i,:);
            H = [H;H_temp];
            %         end
        
        if opt_dist.FLAGS.compare_with_CI
            
            H_temp_CI = zeros(opt_dist.dimObs,size(opt_dist.x_gt,1));
               h_temp_CI = opt_dist.C(idx_agent,:) * x_pred_CI;
                %         if norm(z_temp) <=opt_dist.obs.Range
                id_vis_CI = [id_vis_CI,i];
                h_CI = [h_CI;h_temp_CI];
                H_temp_CI = opt_dist.C(i,:);
                H_CI = [H_CI;H_temp_CI];
                %         end
        end
    end
    opt_dist.result.obs.h{idx_agent} = h;
    opt_dist.result.obs.H{idx_agent} = H;
    opt_dist.result.obs.id_vis{idx_agent} = id_vis;
    
    if opt_dist.FLAGS.compare_with_CI
        opt_dist.result.obs.h_CI{idx_agent} = h_CI;
        opt_dist.result.obs.H_CI{idx_agent} = H_CI;
        opt_dist.result.obs.id_vis_CI{idx_agent} = id_vis_CI;
        
    end
end
end
