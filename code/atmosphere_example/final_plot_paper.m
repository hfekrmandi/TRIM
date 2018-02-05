close all

for i_file=1:11
    clear bcs_perf error_results perf_index range_prob x0
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
        e_bar_=[e_bar_;mean(mean(error_results{i}.e_bar_))];
        e_bar_vs_cen=[e_bar_vs_cen,mean(mean(error_results{i}.e_bar_vs_cen))];
        e_P_bar_cen=[e_P_bar_cen;mean(mean(error_results{i}.e_P_bar_cen))];
        e_P_bar_cen_2=[e_P_bar_cen_2;mean(mean(error_results{i}.e_P_bar_cen_2))];
        e_P_bar_root =[e_P_bar_root;mean(mean(error_results{i}.e_P_bar_root))]
        e_P_CI_root =[e_P_CI_root;mean(mean(error_results{i}.e_P_CI_root))]
        
        
        e_bar_CI=[e_bar_CI;mean(mean(error_results{i}.e_bar_CI))];
        e_bar_CI_vs_cen=[e_bar_CI_vs_cen;mean(mean(error_results{i}.e_bar_CI_vs_cen))];
        e_P_bar_CI_cen=[e_P_bar_CI_cen;mean(mean(error_results{i}.e_P_bar_CI_cen))];
        e_P_bar_CI_cen_2=[e_P_bar_CI_cen_2;mean(mean(error_results{i}.e_P_bar_CI_cen_2))];
        
        e_bar_BC_dist = [e_bar_BC_dist;mean(mean(error_results{i}.e_bar_BC_dist))];
        e_bar_CI_BC_dist = [e_bar_CI_BC_dist;mean(mean(error_results{i}.e_bar_CI_BC_dist))];
        
        con_perc_ci = [con_perc_ci; mean(mean(error_results{i}.con_perc_ci))];
        con_perc_bar = [con_perc_bar; mean(mean(error_results{i}.con_perc_bar))];
        con_perc_cen = [con_perc_cen; mean(mean(error_results{i}.con_perc_cen))];
        
    end
    h_err1 = figure
    h_err2 = figure
    h_err3 = figure
    
    % for i_agent = 1:9
    figure(h_err1)
    subplot(321)
    plot(e_bar_); hold on
    title('e_{bar} = ||x_{bar} - x_{gt}||')
    grid on; grid minor;
    
    subplot(322)
    plot(e_bar_vs_cen); hold on
    title('e_{bar vs. cen}=||x_{bar} - x_{cen}||')
    grid on; grid minor;
    
    subplot(323)
    plot(e_P_bar_cen); hold on
    title('e_{P_{bar} vs. P_{cen}}= det(P_{cen})/det(P_{bar})')
    grid on; grid minor;
    
    subplot(324)
    plot(e_bar_CI); hold on
    title('e_{CI vs. gt}=||x_{CI} - x_{gt}||')
    grid on; grid minor;
    
    subplot(325)
    plot(e_bar_CI_vs_cen); hold on
    title('e_{CI vs. cen}=||x_{CI} - x_{cen}||')
    grid on; grid minor;
    
    subplot(326)
    plot(e_P_bar_CI_cen); hold on
    title('e_{P_{CI} vs. P_{cen}}= det(P_{cen})/det(P_{CI})')
    grid on; grid minor;
    
    figure(h_err2)
    subplot(221)
    plot(e_bar_BC_dist); hold on
    title('BCdist_{bar}')
    grid on; grid minor;
    subplot(222)
    plot(e_bar_CI_BC_dist); hold on
    title('BCdist_{CI}')
    grid on; grid minor;
    subplot(2,2,[3 4])
    plot(e_bar_BC_dist); hold on
    plot(e_bar_CI_BC_dist); hold on
    grid on; grid minor;
    
    
    
    figure(h_err3)
    subplot(221)
    plot(e_P_CI_root); hold on
    title('trace(P_{cen})/trace(P_{bar})')
    grid on; grid minor;
    subplot(222)
    plot(e_P_bar_root); hold on
    title('trace(P_{cen})/trace(P_{CI})')
    grid on; grid minor;
    subplot(2,2,[3 4])
    plot(e_P_bar_cen_2); hold on
    plot(e_P_bar_CI_cen_2); hold on
    grid on; grid minor;
    
    
    % end
    close all
    figure
    subplot(221)
    plot(e_bar_BC_dist); hold on
    plot(e_bar_CI_BC_dist); hold on
    grid on; grid minor;
    legend('Hybrid','CI')
    title('Bhatacharia Distance')
    bc_dist(:,i_file) = [mean(e_bar_BC_dist);mean(e_bar_CI_BC_dist)];
    
    subplot(222)
    plot(e_P_bar_cen_2); hold on
    plot(e_P_bar_CI_cen_2); hold on
    legend('Hybrid','CI')
    title('trace ratio')
    trace_(:,i_file) = [mean(e_P_bar_cen_2);mean(e_P_bar_CI_cen_2)];
    
    subplot(223)
    plot((e_P_bar_cen).^(1/80)); hold on
    plot(e_P_bar_CI_cen.^(1/80)); hold on
    legend('Hybrid','CI')
    title('det ratio')
    det_ratio(:,i_file) = [mean((e_P_bar_cen).^(1/80));mean(e_P_bar_CI_cen.^(1/80))];
    
    subplot(224)
    plot(e_bar_); hold on
    plot(e_bar_CI); hold on
    plot(e_cen_); hold on;
    legend('Hybrid','CI','Centralized')
    title('rmse error w.r.t ground truth')
    rmse_(:,i_file) = [mean(e_bar_);mean(e_bar_CI);mean(e_cen_)];
    
    
    ci(i_file) = mean(mean(con_perc_ci(:)));
    bar_(i_file) = mean(mean(con_perc_bar(:)));
    cen_(i_file)= mean(mean(con_perc_cen(:)));
end
    close all

figure
% prob_range = [0 0.2 0.4 0.5 0.8 1];
prob_range = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

subplot(221)
plot(prob_range,bc_dist(1,:)); hold on;
plot(prob_range,bc_dist(2,:)); hold on;
plot([0 0.1 0.2 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],[1 1 1 0 0 0 0 0 0 0 0 0]); hold on;

    legend('Hybrid','CI','Centralized (wheni network is connected)')
    title('Bhatacharia Distance')
    grid on


subplot(222)
plot(prob_range,trace_(1,:)); hold on;
plot(prob_range,trace_(2,:)); hold on;
    legend('Hybrid','CI')
    title('trace ratio')
    grid on


subplot(223)
plot(prob_range,det_ratio(1,:)); hold on;
plot(prob_range,det_ratio(2,:)); hold on;
    legend('Hybrid','CI')
    title('det ratio')
    grid on

subplot(224)
plot(prob_range,rmse_(1,:)); hold on;
plot(prob_range,rmse_(2,:)); hold on;
plot(prob_range,rmse_(3,:)); hold on;

    legend('Hybrid','CI','Centralized')
    title('rmse error w.r.t ground truth')
    grid on
figure
t = uitable('Data', [ci;bar_;cen_], 'ColumnName', {'0', '0.2', '0.4','0.6','0.8','1'},'RowName', {'CI', 'ConCint', 'Centralized'})
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);   

figure
subplot(211)
% ylim([0.6 1])

plot(prob_range,det_ratio(1,:),'-ob','LineWidth',2); hold on;
plot(prob_range,det_ratio(2,:),'-.sr','LineWidth',2); hold on;
% plot([0 0.1 0.2 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],[1 1 1 0 0 0 0 0 0 0 0 0],'-.g','LineWidth',3); hold on;
ax = gca;

ylim([0.6 1.05])
    grid on
ax.YTick = [0.6:0.1:1];
ax.XTick = [0:0.1:1];
    legend('Hybrid','CI')
    title('{[det ratio]}^{1/80}')
    
    xlabel('Probablity of Link Failure')
    ylabel('{[det ratio]}^{1/80}')
    
subplot(212)
plot(prob_range,bc_dist(1,:),'-ob','LineWidth',2); hold on;
plot(prob_range,bc_dist(2,:),'-.sr','LineWidth',2); hold on;
plot([0 0.1 0.2 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],[1 1 1 0 0 0 0 0 0 0 0 0],'-.g','LineWidth',3); hold on;

ylim([0 1.05])
ax = gca;
ax.YTick = [0:0.2:1];
ax.XTick = [0:0.1:1];

    legend('Hybrid','CI','Centralized')
    title('Bhattacharyya Distance')
    xlabel('Probablity of Link Failure')

    ylabel('D_{B}')
    grid on

