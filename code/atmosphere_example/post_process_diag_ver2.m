function [perf_index,error_results,opt_dist_result,bcs_perf] =post_process_diag_ver2()
global opt_dist
max_it = opt_dist.nSteps;
start_step =1;
for i_consensus=start_step:max_it;
    est{i_consensus - start_step +1}.P_cen = inv(opt_dist.result.est{1}.Y_cen);
    est{i_consensus - start_step +1}.x_cen = est{i_consensus - start_step +1}.P_cen*(opt_dist.result.est{1}.y_cen);
    
    e_cen(i_consensus - start_step +1) = sqrt(immse(est{i_consensus - start_step +1}.x_cen , opt_dist.result.gt.x_bar));
    con_perc_cen(i_consensus - start_step +1)  = consistency_percentage(est{i_consensus - start_step +1}.x_cen ,opt_dist.result.gt.x_bar, est{i_consensus - start_step +1}.P_cen );
    cen_diag = diag(est{i_consensus - start_step +1}.P_cen);
    e_P_cen_root(i_consensus - start_step +1) = (det(est{i_consensus - start_step +1}.P_cen))^(1/80) ;
    
    for j_agent = 1 : opt_dist.nAgents
        if opt_dist.FLAGS.our_method
            est{i_consensus - start_step +1}.P_bar(:,:,j_agent) = inv(opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent));
            est{i_consensus - start_step +1}.x_bar(:,j_agent) = est{i_consensus - start_step +1}.P_bar(:,:,j_agent) *(opt_dist.result.est{i_consensus}.y_bar(:,j_agent));
            e_bar_(i_consensus - start_step +1,j_agent) = sqrt(immse(est{i_consensus - start_step +1}.x_bar(:,j_agent) , opt_dist.result.gt.x_bar));
            e_bar_vs_cen(i_consensus - start_step +1,j_agent) = sqrt(immse(est{i_consensus - start_step +1}.x_bar(:,j_agent) , est{i_consensus - start_step +1}.x_cen ));
            %         e_P_bar_cen(i_consensus - start_step +1,j_agent) = norm(est{i_consensus - start_step +1}.P_bar(:,:,j_agent) - est{i_consensus - start_step +1}.P_cen );
            con_perc_bar(i_consensus - start_step +1,j_agent)  = consistency_percentage(est{i_consensus - start_step +1}.x_bar(:,j_agent),opt_dist.result.gt.x_bar, est{i_consensus - start_step +1}.P_bar(:,:,j_agent) );
%             bar_diag(i_consensus - start_step +1,j_agent,:) = diag( est{i_consensus - start_step +1}.P_bar(:,:,j_agent) );
            
            e_P_bar_root(i_consensus - start_step +1,j_agent) = (det(est{i_consensus - start_step +1}.P_bar(:,:,j_agent)))^(1/80) ;
            
            e_P_bar_cen(i_consensus - start_step +1,j_agent) = det(est{i_consensus - start_step +1}.P_cen)/det(est{i_consensus - start_step +1}.P_bar(:,:,j_agent)) ;
            e_P_bar_cen_2(i_consensus - start_step +1,j_agent) = trace(est{i_consensus - start_step +1}.P_cen)/trace(est{i_consensus - start_step +1}.P_bar(:,:,j_agent)) ;
            
            e_bar_BC_dist(i_consensus - start_step +1,j_agent) = BC_distance(est{i_consensus - start_step +1}.x_cen,est{i_consensus - start_step +1}.P_cen,...
                est{i_consensus - start_step +1}.x_bar(:,j_agent),est{i_consensus - start_step +1}.P_bar(:,:,j_agent));
        end
        if opt_dist.FLAGS.compare_with_CI
            est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent) = inv(opt_dist.result.est{i_consensus}.Y_bar_CI(:,:,j_agent));
%                         bar_CI(i_consensus,j_agent,:) = diag( est{i_consensus}.P_bar_CI(:,:,j_agent) );

            
            est{i_consensus - start_step +1}.x_bar_CI(:,j_agent) = est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent) *(opt_dist.result.est{i_consensus}.y_bar_CI(:,j_agent));
            con_perc_ci(i_consensus - start_step +1,j_agent)  = consistency_percentage(est{i_consensus - start_step +1}.x_bar_CI(:,j_agent),opt_dist.result.gt.x_bar, est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent) );
            e_bar_CI(i_consensus - start_step +1,j_agent) = sqrt(immse(est{i_consensus - start_step +1}.x_bar_CI(:,j_agent) , opt_dist.result.gt.x_bar));
            e_bar_CI_vs_cen(i_consensus - start_step +1,j_agent) = sqrt(immse(est{i_consensus - start_step +1}.x_bar_CI(:,j_agent) , est{i_consensus - start_step +1}.x_cen ));
            e_P_bar_CI_cen(i_consensus - start_step +1,j_agent) = det(est{i_consensus - start_step +1}.P_cen)/det(est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent)) ;
            e_P_bar_CI_cen_2(i_consensus - start_step +1,j_agent) = trace(est{i_consensus - start_step +1}.P_cen)/trace(est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent)) ;
            e_P_CI_root(i_consensus - start_step +1,j_agent) = (det(est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent)))^(1/80) ;
            e_bar_CI_BC_dist(i_consensus - start_step +1,j_agent) = BC_distance(est{i_consensus - start_step +1}.x_cen,est{i_consensus - start_step +1}.P_cen,...
                est{i_consensus - start_step +1}.x_bar_CI(:,j_agent),est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent));
            % e_P_bar_CI_cen(i_consensus - start_step +1,j_agent) = norm(est{i_consensus - start_step +1}.P_bar_CI(:,:,j_agent) - est{i_consensus - start_step +1}.P_cen );
        end
    end
    
    
    
end
error_results.e_cen_  = e_cen;
error_results.e_P_cen_root = e_P_cen_root;

if opt_dist.FLAGS.our_method
    error_results.e_bar_ = e_bar_;
%     error_results.bar_diag = bar_diag;
    error_results.e_bar_vs_cen = e_bar_vs_cen;
    error_results.e_P_bar_cen = e_P_bar_cen;
    error_results.e_P_bar_cen_2 = e_P_bar_cen_2;
    error_results.e_P_bar_root = e_P_bar_root;
    
    error_results.e_bar_BC_dist = e_bar_BC_dist;
end
if opt_dist.FLAGS.compare_with_CI
    error_results.e_bar_CI_BC_dist = e_bar_CI_BC_dist;
%         error_results.bar_CI = bar_CI;

    error_results.e_bar_CI = e_bar_CI;
    error_results.e_bar_CI_vs_cen = e_bar_CI_vs_cen;
    error_results.e_P_bar_CI_cen = e_P_bar_CI_cen;
    error_results.e_P_bar_CI_cen_2 = e_P_bar_CI_cen_2;
    
    error_results.e_P_CI_root = e_P_CI_root;
    
    
    
    error_results.con_perc_ci = con_perc_ci;
end
error_results.con_perc_bar = con_perc_bar;
error_results.con_perc_cen = con_perc_cen;



opt_dist.result.prior.x_cen = est{i_consensus - start_step +1}.x_cen ;
opt_dist.result.prior.P_cen = est{i_consensus - start_step +1}.P_cen ;

for i_agent = 1 : opt_dist.nAgents
    if opt_dist.FLAGS.compare_with_CI
        opt_dist.result.prior.x_bar_CI(:,i_agent) = est{end}.x_bar_CI(:,i_agent);
        opt_dist.result.prior.P_bar_CI(:,:,i_agent) = est{end}.P_bar_CI(:,:,i_agent);
    end
    if opt_dist.FLAGS.our_method
        
        opt_dist.result.prior.x_bar(:,i_agent) = est{end}.x_bar(:,i_agent);
        opt_dist.result.prior.P_bar(:,:,i_agent) = est{end}.P_bar(:,:,i_agent);
    end
end
opt_dist.dataG.Graph
opt_dist.result.graph = opt_dist.dataG.Graph;
% est_ = [est{[1:5,65:opt_dist.nSteps ]}];
% xgt= opt_dist.result.gt.x_bar;

% perf_index = [mean(e_P_bar_cen(end,:)),mean(e_P_bar_CI_cen(end,:))];
perf_index = [];
bcs_perf = [];
if opt_dist.FLAGS.our_method
    
    perf_index = [perf_index,mean(e_P_bar_cen(:))];
    bcs_perf = [bcs_perf,mean(error_results.e_bar_BC_dist(:))];
end
if opt_dist.FLAGS.compare_with_CI
    
    perf_index = [perf_index,mean(e_P_bar_CI_cen(:))];
        bcs_perf = [bcs_perf,mean(error_results.e_bar_CI_BC_dist(:))];

end
if opt_dist.FLAGS.debug
    opt_dist.figures.fig_cov_debug
end
if any(opt_dist.iter_interest==opt_dist.i_step)
    opt_dist_result = est;
else
    opt_dist_result = [];
end



