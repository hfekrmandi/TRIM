function consenus()
global opt_dist

% initialize the concsensus values
for i_agent = 1 : opt_dist.nAgents
    if opt_dist.FLAGS.our_method
        opt_dist.result.consenus{i_agent}.Y_prior{1} = opt_dist.result.pred.Y_bar(:,:,i_agent) ;
        opt_dist.result.consenus{i_agent}.y_prior{1} = opt_dist.result.pred.y_bar(:,i_agent);
        H = opt_dist.result.obs.H{i_agent};
        z = opt_dist.sim.obs.z{i_agent};
        opt_dist.result.consenus{i_agent}.delta_I{1} = 1/(opt_dist.sim.obs.r_var{i_agent})*(H'*H);
        opt_dist.result.consenus{i_agent}.delta_i{1} =  1/(opt_dist.sim.obs.r_var{i_agent})*H'*z;
        opt_dist.result.consenus{i_agent}.group_set_ci{1} = i_agent;

    end
    if opt_dist.FLAGS.compare_with_CI
        
        opt_dist.result.consenus{i_agent}.Y_prior_CI{1} = opt_dist.result.pred.Y_bar_CI(:,:,i_agent) ;
        opt_dist.result.consenus{i_agent}.y_prior_CI{1} = opt_dist.result.pred.y_bar_CI(:,i_agent);
        H_CI = opt_dist.result.obs.H_CI{i_agent};
        z_CI = opt_dist.sim.obs.z_CI{i_agent};
        opt_dist.result.consenus{i_agent}.delta_I_CI{1} = 1/(opt_dist.sim.obs.r_CI_var{i_agent})*(H_CI'*H_CI) +opt_dist.result.pred.Y_bar_CI(:,:,i_agent);
        opt_dist.result.consenus{i_agent}.delta_i_CI{1} =  1/(opt_dist.sim.obs.r_CI_var{i_agent})*H_CI'*z_CI+ opt_dist.result.pred.y_bar_CI(:,i_agent);
        opt_dist.result.consenus{i_agent}.group_set{1} = i_agent;
    end
    
end


% [av_delta_I,av_delta_i] = calculate_aver_info(opt_dist.result);

%% centralized update
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

opt_dist.result.est{1}.Y_cen = opt_dist.result.pred.Y_cen + opt_dist.result.update.delta_I_cen;
opt_dist.result.est{1}.y_cen = opt_dist.result.pred.y_cen + opt_dist.result.update.delta_i_cen;
% for s_agent = 1 : opt_dist.nAgents
%     for s_agent = 1 : opt_dist.nAgents
%
%     end
% end

i_consensus = 1;
%% our method first step update
if opt_dist.FLAGS.our_method
    for j_agent = 1 : opt_dist.nAgents
        
        opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent) = opt_dist.result.consenus{j_agent}.Y_prior{i_consensus}+...
           numel(opt_dist.result.consenus{i_agent}.group_set{1})*opt_dist.result.consenus{j_agent}.delta_I{i_consensus};
        opt_dist.result.est{i_consensus}.y_bar(:,j_agent) = opt_dist.result.consenus{j_agent}.y_prior{i_consensus}+...
            numel(opt_dist.result.consenus{i_agent}.group_set{1})*opt_dist.result.consenus{j_agent}.delta_i{i_consensus};
    end
end
if  opt_dist.FLAGS.compare_with_CI
    for j_agent = 1 : opt_dist.nAgents
        
        opt_dist.result.est{i_consensus}.Y_bar_CI(:,:,j_agent) = opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus};
        opt_dist.result.est{i_consensus}.y_bar_CI(:,j_agent) = opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus};
    end
end



% i_consensus = 1;
% is_converged = 0;
% while ~is_converged
if opt_dist.FLAGS.verbose
    h_wait = waitbar(0,'Consensus process starts..');
end
for i_consensus = 2:opt_dist.nSteps
    %  i_consensus = i_consensus + 1;
    opt_dist.result.est{i_consensus}.Y_cen = opt_dist.result.pred.Y_cen + opt_dist.result.update.delta_I_cen;
    opt_dist.result.est{i_consensus}.y_cen = opt_dist.result.pred.y_cen + opt_dist.result.update.delta_i_cen;
    for j_agent = 1 : opt_dist.nAgents
        if opt_dist.FLAGS.our_method
            opt_dist.result.consenus{j_agent}.Y_prior{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.Y_prior{i_consensus-1}));
            opt_dist.result.consenus{j_agent}.y_prior{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.y_prior{i_consensus-1}));
            opt_dist.result.consenus{j_agent}.delta_I{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.delta_I{i_consensus-1}));
            opt_dist.result.consenus{j_agent}.delta_i{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.delta_i{i_consensus-1}));
            I_local=[];P_local=[];i_local = [];prev_delta_I=[];prev_delta_i=[];
            
        end
        if opt_dist.FLAGS.compare_with_CI
            opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus-1}));
            opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus-1}));
            opt_dist.result.consenus{j_agent}.delta_I_CI{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.delta_I_CI{i_consensus-1}));
            opt_dist.result.consenus{j_agent}.delta_i_CI{i_consensus} = zeros(size(opt_dist.result.consenus{j_agent}.delta_i_CI{i_consensus-1}));
            I_local_CI=[];P_local_CI=[];i_local_CI=[];
            
        end
        
        %% prior weight calc
        idx_neighbours = find(opt_dist.dataG.Graph.Adj(j_agent,:));
        
        for j_neigh=1:length(idx_neighbours)
            if opt_dist.FLAGS.our_method
                I_local(:,:,j_neigh) = (opt_dist.result.consenus{idx_neighbours(j_neigh)}.Y_prior{i_consensus-1});
                i_local(:,:,j_neigh) = (opt_dist.result.consenus{idx_neighbours(j_neigh)}.y_prior{i_consensus-1});
                P_local(:,:,j_neigh) = inv(I_local(:,:,j_neigh));
                prev_delta_I(:,:,j_neigh) = opt_dist.result.consenus{j_neigh}.delta_I{i_consensus-1};
                prev_delta_i(:,:,j_neigh) = opt_dist.result.consenus{j_neigh}.delta_i{i_consensus-1};
            end
            if opt_dist.FLAGS.compare_with_CI
                I_local_CI(:,:,j_neigh) = opt_dist.result.consenus{idx_neighbours(j_neigh)}.delta_I_CI{i_consensus-1};
                i_local_CI(:,:,j_neigh) = (opt_dist.result.consenus{idx_neighbours(j_neigh)}.delta_i_CI{i_consensus-1});
                P_local_CI(:,:,j_neigh) = inv(I_local_CI(:,:,j_neigh));
            end
        end
        if opt_dist.FLAGS.our_method
            %             weights_ci =calc_ci_weights(P_local,'det');
            log_message('our method CI')
            [weights_ci,inf_mat,inf_vect] =calc_ci_weights_ver2(P_local,i_local,'det');
            opt_dist.result.consenus{j_agent}.Y_prior{i_consensus} = inf_mat;
            opt_dist.result.consenus{j_agent}.y_prior{i_consensus} = inf_vect;
            %             mh_weights_row = opt_dist.dataG.Graph.p(j_agent,idx_neighbours);
            %             opt_dist.result.consenus{j_agent}.delta_I{i_consensus} = special_dot_sum(mh_weights_row,prev_delta_I,0);
            %             opt_dist.result.consenus{j_agent}.delta_i{i_consensus} =special_dot_sum(mh_weights_row,prev_delta_i,0);
            opt_dist.result.consenus{j_agent}.group_set{i_consensus}  = find( opt_dist.dataG.Graph.p(j_agent,:));
            updated_delta_I = []; updated_delta_i=[];
            for k_agent = 1 : opt_dist.nAgents
                p_jk = opt_dist.dataG.Graph.p(j_agent,k_agent);
                
                %             b_jk = covariance_intersection()
                updated_delta_I = opt_dist.result.consenus{j_agent}.delta_I{i_consensus} + p_jk*opt_dist.result.consenus{k_agent}.delta_I{i_consensus-1};
                updated_delta_i = opt_dist.result.consenus{j_agent}.delta_i{i_consensus} + p_jk*opt_dist.result.consenus{k_agent}.delta_i{i_consensus-1};
                opt_dist.result.consenus{j_agent}.delta_I{i_consensus} = updated_delta_I;
                opt_dist.result.consenus{j_agent}.delta_i{i_consensus} =updated_delta_i;
                opt_dist.result.consenus{j_agent}.group_set{i_consensus} = union( opt_dist.result.consenus{j_agent}.group_set{i_consensus},...
                                                                               opt_dist.result.consenus{k_agent}.group_set{i_consensus-1});
            end
            if opt_dist.dataG.is_connected
                diff_I = max(max((av_delta_I - opt_dist.result.consenus{j_agent}.delta_I{i_consensus})))
                diff_i = max((av_delta_i - opt_dist.result.consenus{j_agent}.delta_i{i_consensus}))
            end
            opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent) = opt_dist.result.consenus{j_agent}.Y_prior{i_consensus}+...
                numel(opt_dist.result.consenus{j_agent}.group_set{i_consensus})*opt_dist.result.consenus{j_agent}.delta_I{i_consensus};
            opt_dist.result.est{i_consensus}.y_bar(:,j_agent) = opt_dist.result.consenus{j_agent}.y_prior{i_consensus}+...
                 numel(opt_dist.result.consenus{j_agent}.group_set{i_consensus})*opt_dist.result.consenus{j_agent}.delta_i{i_consensus};
            
        end
        if opt_dist.FLAGS.compare_with_CI
            log_message('Pure CI')
            [weights_ci_all,inf_mat_ci,inf_vect_ci] =calc_ci_weights_ver2(P_local_CI,(i_local_CI),'det');
            opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus} = inf_mat_ci;
            opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus} = inf_vect_ci;
            opt_dist.result.consenus{j_agent}.delta_I_CI{i_consensus} = inf_mat_ci;
            opt_dist.result.consenus{j_agent}.delta_i_CI{i_consensus} = inf_vect_ci;
            opt_dist.result.est{i_consensus}.Y_bar_CI(:,:,j_agent) = opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus};
            opt_dist.result.est{i_consensus}.y_bar_CI(:,j_agent) = opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus};
        end
    end
    if opt_dist.FLAGS.verbose
        waitbar(i_consensus / opt_dist.nSteps, h_wait,['step = ',num2str(opt_dist.i_step),' consensus iteration = ',num2str(i_consensus)])
    end
end
if opt_dist.FLAGS.verbose
    
    close(h_wait)
end
% if opt_dist.FLAGS.debug
%     opt_dist.figures.fig_cov_debug
% end

end


%% previous code
%         for k_priori=1:size(weights_ci,1)
%             if opt_dist.FLAGS.our_method
%                 b_jk = weights_ci(k_priori);
%                 updated_Y_prior = opt_dist.result.consenus{j_agent}.Y_prior{i_consensus} + b_jk*opt_dist.result.consenus{idx_neighbours(k_priori)}.Y_prior{i_consensus-1};
%                 updated_y_prior = opt_dist.result.consenus{j_agent}.y_prior{i_consensus} + b_jk*opt_dist.result.consenus{idx_neighbours(k_priori)}.y_prior{i_consensus-1};
%                 opt_dist.result.consenus{j_agent}.Y_prior{i_consensus} = updated_Y_prior;
%                 opt_dist.result.consenus{j_agent}.y_prior{i_consensus} = updated_y_prior;
%             end
%             if opt_dist.FLAGS.compare_with_CI
%                     b_jk_all = weights_ci_all(k_priori);
%                 updated_delta_I_CI = opt_dist.result.consenus{j_agent}.delta_I_CI{i_consensus} + b_jk_all*I_local_CI(:,:,k_priori) ;
%                 updated_delta_i_CI = opt_dist.result.consenus{j_agent}.delta_i_CI{i_consensus} + b_jk_all*i_local_CI(:,:,k_priori);
%                 opt_dist.result.consenus{j_agent}.delta_I_CI{i_consensus} = updated_delta_I_CI;
%                 opt_dist.result.consenus{j_agent}.delta_i_CI{i_consensus} = updated_delta_i_CI;
%
%                 opt_dist.result.consenus{j_agent}.Y_prior_CI{i_consensus} = updated_delta_I_CI;
%                 opt_dist.result.consenus{j_agent}.y_prior_CI{i_consensus} = updated_delta_i_CI;
%             end
%         end
%         if opt_dist.FLAGS.our_method
%             for k_agent = 1 : opt_dist.nAgents
%                 p_jk = opt_dist.dataG.Graph.p(j_agent,k_agent);
%
%                 %             b_jk = covariance_intersection()
%                 updated_delta_I = opt_dist.result.consenus{j_agent}.delta_I{i_consensus} + p_jk*opt_dist.result.consenus{k_agent}.delta_I{i_consensus-1};
%                 updated_delta_i = opt_dist.result.consenus{j_agent}.delta_i{i_consensus} + p_jk*opt_dist.result.consenus{k_agent}.delta_i{i_consensus-1};
%                 opt_dist.result.consenus{j_agent}.delta_I{i_consensus} = updated_delta_I;
%                 opt_dist.result.consenus{j_agent}.delta_i{i_consensus} =updated_delta_i;
%
%             end
%         end