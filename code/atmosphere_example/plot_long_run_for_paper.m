clear all
close all
n_dim = 150;
for i_file=1:11
    clear bcs_perf error_results perf_index range_prob x0
%     file_name = ['atmosh_reduced_order',num2str(i_file)];
    file_name = ['file__0',num2str(i_file)];

    load(file_name)
    
    e_bar_=[];
    e_cen_=[];
    
    e_bar_vs_cen=[];
    e_P_bar_cen=[];
    e_bar_CI=[];
    e_bar_CI_vs_cen=[];
    e_P_bar_CI_cen=[];
    e_bar_BC_dist = [];
    e_bar_CI_BC_dist = [];
    con_perc_ci = [];
    con_perc_bar = [];
    con_perc_cen = [];
    e_P_bar_cen_2 = [];
    e_P_bar_CI_cen_2 = [];
    e_P_bar_root =[]
    e_P_CI_root =[];
    
    for i=1:numel(error_results)
        e_cen_ = [e_cen_,mean(error_results{i}.e_cen_)];
        e_bar_=[e_bar_;mean(mean(error_results{i}.e_bar_,2))];
        e_bar_vs_cen=[e_bar_vs_cen,mean(mean(error_results{i}.e_bar_vs_cen,2))];
        e_P_bar_cen=[e_P_bar_cen;mean(mean(error_results{i}.e_P_bar_cen,2))];
        e_P_bar_cen_2=[e_P_bar_cen_2;mean(mean(error_results{i}.e_P_bar_cen_2,2))];
        e_P_bar_root =[e_P_bar_root;mean(mean(error_results{i}.e_P_bar_root,2))]
        e_P_CI_root =[e_P_CI_root;mean(mean(error_results{i}.e_P_CI_root,2))]
        
        
        e_bar_CI=[e_bar_CI;mean(mean(error_results{i}.e_bar_CI,2))];
        e_bar_CI_vs_cen=[e_bar_CI_vs_cen;mean(mean(error_results{i}.e_bar_CI_vs_cen,2))];
        e_P_bar_CI_cen=[e_P_bar_CI_cen;mean(mean(error_results{i}.e_P_bar_CI_cen,2))];
        e_P_bar_CI_cen_2=[e_P_bar_CI_cen_2;mean(mean(error_results{i}.e_P_bar_CI_cen_2,2))];
        
        e_bar_BC_dist = [e_bar_BC_dist;mean(mean(error_results{i}.e_bar_BC_dist,2))];
        e_bar_CI_BC_dist = [e_bar_CI_BC_dist;mean(mean(error_results{i}.e_bar_CI_BC_dist,2))];
        
        con_perc_ci = [con_perc_ci; mean(mean(error_results{i}.con_perc_ci))];
        con_perc_bar = [con_perc_bar; mean(mean(error_results{i}.con_perc_bar))];
        con_perc_cen = [con_perc_cen; mean(mean(error_results{i}.con_perc_cen))];
        
    end
    
    
    
 
    close all
    figure
    subplot(313)
    plot(e_bar_BC_dist,'-r','LineWidth',2); hold on
    plot(e_bar_CI_BC_dist,'--b','LineWidth',2); hold on
    grid on; grid minor;
    xlabel('steps')
   ylabel('D_{B}')

    legend('Hybrid','CI')
    title('Bhatacharia Distance')
        bc_dist(:,i_file) = [mean(e_bar_BC_dist);mean(e_bar_CI_BC_dist)];
    
    
    subplot(312)
    plot((e_P_bar_cen).^(1/n_dim),'-r','LineWidth',2); hold on
    plot(e_P_bar_CI_cen.^(1/n_dim),'--b','LineWidth',2); hold on
    legend('Hybrid','CI')
    xlabel('steps')
       ylabel('{[det ratio]^{1/n_dim}}')

    title('{[det ratio]^{1/n_dim}}')
        grid on; grid minor;

        det_ratio(:,i_file) = [mean((e_P_bar_cen).^(1/n_dim));mean(e_P_bar_CI_cen.^(1/n_dim))];
    
    subplot(311)
    plot(e_bar_,'-r','LineWidth',2); hold on
    plot(e_bar_CI,'--b','LineWidth',2); hold on
    plot(e_cen_,'.-g','LineWidth',2); hold on;
    xlabel('steps')
    
    legend('Hybrid','CI','Centralized')
    ylabel('rmse error')
    title('rmse error')
        grid on; grid minor;

%         rmse_(:,i_file) = [mean(e_bar_);mean(e_bar_CI);(e_cen_')];
end 
close all
figure
  plot(bc_dist(1,:),'LineWidth',2,'Marker','o','MarkerSize',7); hold on;
  plot(bc_dist(2,:),'LineWidth',2,'Marker','+','MarkerSize',7); hold on;
% plot(bc_dist(1,:),'LineWidth',2); hold on;
%   plot(bc_dist(2,:),'LineWidth',2); hold on;
%   
legend('Hybrid','CI')
  
  xlabel('Probability of Link Failure','FontSize',16)
  ylabel('D_B','FontSize',16)
title('Bhattacharyya distance','FontSize',16)    
ax = gca;
ax.XTick = [0 1 2 3 4 5 6];
% ax.YTick = [-1 -0.5 0 0.5 1];

ax.XTickLabel = {'0.0','0.2','0.4','0.6','0.8','1'};
% ax.YTickLabel = {'min = -1','-0.5','0','0.5','max = 1'};
set(ax,'FontSize',14)
grid on
grid minor