function post_process()
global opt_dist
persistent error_results
max_it = 70;
% opt_dist.result.est.P_cen = inv(opt_dist.result.est{max_it}.Y_cen );
% opt_dist.result.est.x_cen = opt_dist.result.est.P_cen*opt_dist.result.est.y_cen ;

% opt_dist.result.est.e_cen = opt_dist.result.est.x_cen - opt_dist.result.gt.x_bar;
% opt_dist.result.est.e_cen_norm = norm(opt_dist.result.est.e_cen);

% decentralized update
% for i_consensus=1:max_it;
%     for j_agent = 1 : opt_dist.nAgents
%         opt_dist.result.est{i_consensus}.P_bar(:,:,j_agent) = inv( opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent));
%         opt_dist.result.est{i_consensus}.x_bar(:,j_agent) =   opt_dist.result.est{i_consensus}.P_bar(:,:,j_agent)  * opt_dist.result.est{i_consensus}.y_bar(:,:,j_agent);
%         opt_dist.result.est{i_consensus}.e_bar(:,j_agent) = opt_dist.result.est{i_consensus}.x_bar(:,j_agent) - opt_dist.result.est.x_cen ;
%         opt_dist.result.est{i_consensus}.e_bar_norm(j_agent) = norm(opt_dist.result.est{i_consensus}.e_bar(:,j_agent) );
%         opt_dist.result.est{i_consensus}.E_bar(:,:,j_agent) = opt_dist.result.est{i_consensus}.P_bar(:,j_agent) - opt_dist.result.est.P_cen ;
%         opt_dist.result.est{i_consensus}.E_bar_norm(j_agent) = norm(opt_dist.result.est{i_consensus}.E_bar(:,j_agent) );
%     end
% end
% 
% 
% 
% for i_consensus=1:70;
%     for j_agent = 1 : opt_dist.nAgents
% 
%         e_bar_norm(i_consensus,j_agent)
% 
%         opt_dist.result.est{i_consensus}.e_bar_norm(j_agent) = norm(opt_dist.result.est{i_consensus}.e_bar(:,j_agent) );
%         opt_dist.result.est{i_consensus}.E_bar_norm(:,:,j_agent) = norm(opt_dist.result.est{i_consensus}.E_bar(:,j_agent) );
%     end
% end

for i_consensus=1:70;
    est{i_consensus}.P_cen = inv(opt_dist.result.est{i_consensus}.Y_cen);
    est{i_consensus}.x_cen = est{i_consensus}.P_cen*(opt_dist.result.est{i_consensus}.y_cen);
% % %     e_cen(i_consensus) = sqrt(immse(est{i_consensus}.x_cen , opt_dist.result.gt.x_bar));
% % %     
% % %     
% % %     for j_agent = 1 : opt_dist.nAgents
% % %         est{i_consensus}.P_bar(:,:,j_agent) = inv(opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent));
% % %         est{i_consensus}.x_bar(:,j_agent) = est{i_consensus}.P_bar(:,:,j_agent) *(opt_dist.result.est{i_consensus}.y_bar(:,j_agent));
% % %         e_bar_(i_consensus,j_agent) = norm(est{i_consensus}.x_bar(:,j_agent) - opt_dist.result.gt.x_bar);
% % %         e_bar_vs_cen(i_consensus,j_agent) = norm(est{i_consensus}.x_bar(:,j_agent) - est{i_consensus}.x_cen );
% % %         e_P_bar_cen(i_consensus,j_agent) = norm(est{i_consensus}.P_bar(:,:,j_agent) - est{i_consensus}.P_cen );
% % % 
% % %         if opt_dist.FLAGS.compare_with_CI
% % %             est{i_consensus}.P_bar_CI(:,:,j_agent) = inv(opt_dist.result.est{i_consensus}.Y_bar_CI(:,:,j_agent));
% % %             est{i_consensus}.x_bar_CI(:,j_agent) = est{i_consensus}.P_bar_CI(:,:,j_agent) *(opt_dist.result.est{i_consensus}.y_bar_CI(:,j_agent));
% % %             
% % %             e_bar_CI(i_consensus,j_agent) = norm(est{i_consensus}.x_bar_CI(:,j_agent) - opt_dist.result.gt.x_bar);
% % %             e_bar_CI_vs_cen(i_consensus,j_agent) = norm(est{i_consensus}.x_bar_CI(:,j_agent) - est{i_consensus}.x_cen );
% % %             e_P_bar_CI_cen(i_consensus,j_agent) = norm(est{i_consensus}.P_bar_CI(:,:,j_agent) - est{i_consensus}.P_cen );
% % %         end
% % %     end


    e_cen(i_consensus) = sqrt(immse(est{i_consensus}.x_cen , opt_dist.result.gt.x_bar));
    
    
    for j_agent = 1 : opt_dist.nAgents
        est{i_consensus}.P_bar(:,:,j_agent) = inv(opt_dist.result.est{i_consensus}.Y_bar(:,:,j_agent));
        est{i_consensus}.x_bar(:,j_agent) = est{i_consensus}.P_bar(:,:,j_agent) *(opt_dist.result.est{i_consensus}.y_bar(:,j_agent));
        e_bar_(i_consensus,j_agent) = sqrt(immse(est{i_consensus}.x_bar(:,j_agent) , opt_dist.result.gt.x_bar));
        e_bar_vs_cen(i_consensus,j_agent) = sqrt(immse(est{i_consensus}.x_bar(:,j_agent) , est{i_consensus}.x_cen ));
        e_P_bar_cen(i_consensus,j_agent) = norm(est{i_consensus}.P_bar(:,:,j_agent) - est{i_consensus}.P_cen );

        if opt_dist.FLAGS.compare_with_CI
            est{i_consensus}.P_bar_CI(:,:,j_agent) = inv(opt_dist.result.est{i_consensus}.Y_bar_CI(:,:,j_agent));
            est{i_consensus}.x_bar_CI(:,j_agent) = est{i_consensus}.P_bar_CI(:,:,j_agent) *(opt_dist.result.est{i_consensus}.y_bar_CI(:,j_agent));
            
            e_bar_CI(i_consensus,j_agent) = sqrt(immse(est{i_consensus}.x_bar_CI(:,j_agent) , opt_dist.result.gt.x_bar));
            e_bar_CI_vs_cen(i_consensus,j_agent) = sqrt(immse(est{i_consensus}.x_bar_CI(:,j_agent) , est{i_consensus}.x_cen ));
            e_P_bar_CI_cen(i_consensus,j_agent) = norm(est{i_consensus}.P_bar_CI(:,:,j_agent) - est{i_consensus}.P_cen );
        end
    end



end
if isempty(error_results)
    error_results.e_cen = [e_cen'];
    error_results.e_bar = [e_bar_];
    error_results.e_bar_vs_cen = e_bar_vs_cen;
    error_results.e_P_bar_cen = e_P_bar_cen;
    if opt_dist.FLAGS.compare_with_CI
        error_results.e_bar_CI = [e_bar_CI];
        error_results.e_bar_CI_vs_cen = e_bar_CI_vs_cen;
        error_results.e_P_bar_CI_cen = e_P_bar_CI_cen;
    end
else
    error_results.e_cen = [error_results.e_cen;e_cen'];
    error_results.e_bar = [error_results.e_bar;e_bar_];
    error_results.e_bar_vs_cen = [error_results.e_bar_vs_cen;e_bar_vs_cen];
    error_results.e_P_bar_cen = [error_results.e_P_bar_cen;e_P_bar_cen];
    if opt_dist.FLAGS.compare_with_CI
        error_results.e_bar_CI = [error_results.e_bar_CI;e_bar_CI];
        error_results.e_bar_CI_vs_cen = [error_results.e_bar_CI_vs_cen;e_bar_CI_vs_cen];
        error_results.e_P_bar_CI_cen = [error_results.e_P_bar_CI_cen;e_P_bar_CI_cen];
    end
end


close all
agnts_idx = [1 3 5];
figure
subplot(411)
plot(error_results.e_cen) 
title('e_{cen} = norm(x_{cen} - x_{gt})' )

subplot(412)
for j_agent=agnts_idx
plot( error_results.e_bar(:,j_agent) ,'b') ; hold on;
plot( error_results.e_bar_CI(:,j_agent) ,'r') ; hold on;

end
title('e_{bar} =  norm(x_{bar} - x_{gt})')
% legend('1','3','3','4','5')

subplot(413)
for j_agent=agnts_idx
plot( error_results.e_bar_vs_cen(:,j_agent) ,'b') ; hold on;
plot( error_results.e_bar_CI_vs_cen(:,j_agent) ,'r') ; hold on;

end

title('e_{bar\_vs\_cen} = norm(x_{bar} - x_{cen})')
% legend('1','2','3','4','5')

subplot(414)
for j_agent=agnts_idx
plot( error_results.e_P_bar_cen(:,j_agent) ,'b') ; hold on;
plot( error_results.e_P_bar_CI_cen(:,j_agent) ,'r') ; hold on;
end

title('e_{P_{bar\_cen}} = norm(P_{bar} - P_{cen})')
% legend('1','2','3','4','5')



opt_dist.result.prior.x_cen = est{i_consensus}.x_cen ;
opt_dist.result.prior.P_cen = est{i_consensus}.P_cen ;

for i_agent = 1 : opt_dist.nAgents
    if opt_dist.FLAGS.compare_with_CI
        opt_dist.result.prior.x_bar_CI(:,i_agent) = est{max_it}.x_bar_CI(:,i_agent);
        opt_dist.result.prior.P_bar_CI(:,:,i_agent) = est{max_it}.P_bar_CI(:,:,i_agent);
    end
    opt_dist.result.prior.x_bar(:,i_agent) = est{max_it}.x_bar(:,i_agent);
    opt_dist.result.prior.P_bar(:,:,i_agent) = est{max_it}.P_bar(:,:,i_agent);
    
end
est_ = [est{[1:5,65:70]}];
xgt= opt_dist.result.gt.x_bar;
filename_ = sprintf(fullfile(opt_dist.save_dir,'file%d.mat'), opt_dist.i_step);
save(filename_,'est_','xgt')






