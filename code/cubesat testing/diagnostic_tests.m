function diagnostic_tests(error_results,opt_dist_result)
   for i_consensus=1:11
        
        x_cen = est_(i_consensus).x_cen;
        std_cen = sqrt(diag(est_(i_consensus).P_cen));
        e_cen_gt = x_cen-x_gt;
        ow_cen = find(abs(e_cen_gt)>=2*std_cen);
        if ~isempty(ow_cen)
            disp(['step = ',num2str(i_step),'i_consensus = ',num2str(i_consensus),'------------'])
            disp('oo cen')
            disp(ow_cen')
            hist_cen = [hist_cen,ow_cen'];
            disp('-----------')
        end
        for i_agent=1:9
            x_bar = est_(i_consensus).x_bar(:,i_agent);
            std_bar = sqrt(diag(est_(i_consensus).P_bar(:,:,i_agent)));
            
            trace_cov_bar{i_agent} = [trace_cov_bar{i_agent},trace(est_(i_consensus).P_bar(:,:,i_agent) - est_(i_consensus).P_cen)];
            
            x_ci = est_(i_consensus).x_bar_CI(:,i_agent);
            std_ci = sqrt(diag(est_(i_consensus).P_bar_CI(:,:,i_agent)));
            trace_cov_ci{i_agent} = [trace_cov_ci{i_agent}, trace(est_(i_consensus).P_bar_CI(:,:,i_agent)  - est_(i_consensus).P_cen)];
            
            e_ci_gt = x_ci -x_gt;
            e_bar_gt = x_bar - x_gt;
            error_{i_agent} = [error_{i_agent},e_bar_gt];
            std_{i_agent} = [std_{i_agent},std_bar];
            
            
            ow_ci = find(abs(e_ci_gt)>=2*std_ci);
            if ~isempty(ow_ci)
                disp(['step = ',num2str(i_step),'i_consensus = ',num2str(i_consensus),'-------- i_agent = ',num2str(i_agent)'])
                
                disp('oo ci')
                disp(ow_ci')
                disp(e_ci_gt(ow_ci)')
                disp(2*std_ci(ow_ci)')
                disp('-----------')
                hist_ci{i_agent} = [hist_ci{i_agent},ow_ci'];
                
            end
            ow_bar = find(abs(e_bar_gt)>=2*std_bar);
            if ~isempty(ow_bar)
                disp(['step = ',num2str(i_step),'i_consensus = ',num2str(i_consensus),'-------- i_agent = ',num2str(i_agent)'])
                
                disp('oo bar')
                disp(ow_bar')
                
                disp(e_bar_gt(ow_bar)')
                disp(2*std_bar(ow_bar)')
                disp('-----------')
                hist_bar{i_agent} = [hist_bar{i_agent},ow_bar'];
                
            end
            
        end
    end



end