clear all
clc
close all
load('reg_shell2')
figure
prob_ran = [0 0.2 0.4 0.6 0.8 1]
cut_off = [0.01 0.09 0.18 0.32];
reg_ = [2 4 6 8]
for i_reg=1:4
    subplot(2,2,i_reg)
    data_ = squeeze(perf_index_prob(i_reg,:,:));
    h(i_reg) = plot(prob_ran,data_(:,1),'LineWidth',3); hold on
    leg{i_reg} = ['Hybrid, regularity degree = ',num2str(reg_(i_reg) )]
   h_(i_reg) = plot(prob_ran,data_(:,2),'LineWidth',3); hold on
leg_{i_reg} = ['CI, regularity degree = ',num2str(reg_(i_reg) )]


   h__(i_reg) = plot([0 cut_off(i_reg) cut_off(i_reg) 1],[1 1 0 0],'LineWidth',3); hold on
leg__{i_reg} = ['centralized, regularity degree = ',num2str(reg_(i_reg) )]


legend([ h(i_reg) , h_(i_reg),h__(i_reg) ],{leg{i_reg},leg_{i_reg},leg__{i_reg}})
xlabel('Link Fail Prob.')
ylabel('Est.  Performance')
grid on
grid minor
title('Bhattacharyya Distance')
end

figure
for i_reg=1:4
%     subplot(2,2,i_reg)
    data_ = squeeze(perf_index_prob(i_reg,:,:));
    h(i_reg) = plot(prob_ran,data_(:,1),'b'); hold on
    text(0.6,data_(4,1),['\leftarrow',num2str(reg_(i_reg))],'Color','blue','FontSize',12)
%     leg{i_reg} = ['Hybrid, regularity degree = ',num2str(reg_(i_reg) )]
   h_(i_reg) = plot(prob_ran,data_(:,2),'r'); hold on
%        text(0.8,data_(5,2),num2str(reg_(i_reg)),['\leftarrow',num2str(reg_(i_reg))],'Color','red','FontSize',12)
    text(0.8,data_(5,2),['\leftarrow',num2str(reg_(i_reg))],'Color','red','FontSize',12)

% leg_{i_reg} = ['CI, regularity degree = ',num2str(reg_(i_reg) )]
% legend([ h(i_reg) , h_(i_reg) ],{leg{i_reg},leg_{i_reg}})

end
leg_{i_reg} = ['CI'];
leg{i_reg} = ['Hybrid'];
legend([ h(i_reg) , h_(i_reg) ],{leg{i_reg},leg_{i_reg}})

xlabel('Link Fail Prob.')
ylabel('Est.  Performance')
grid on
grid minor
title('Bhattacharyya Distance')




