function agentsVserror_plot()
clearvars
close all
clc

load('cubesat_test_3.mat')
nAgents = opt_dist.nAgents;
agents = 1:nAgents;
fin_step = opt_dist.i_step;
time = 1:fin_step;

T = size(x_gt);
gt_data = zeros(fin_step, nAgents);
%gold_data = zeros(fin_step, nAgents);
hyb_data = zeros(fin_step, nAgents);
ici_data = zeros(fin_step, nAgents);
for i=1:fin_step
    gt_data(i,:) = x_gt{:,i};
    %gold_data(i,:) = x_gold{i,:};
    hyb_data(i,:) = x_Hybrid{i,:};
    ici_data(i,:) = x_ICI{i,:};
end
%%
fig = figure(1);
colors = ['k', 'b', 'r', 'g'];
for i = 1:nAgents
    plot(time, gt_data(:,i)', colors(1));
    hold on
    %plot(time, gold_data(:,i), colors(2));
    plot(time, hyb_data(:,i), colors(3));
    plot(time, ici_data(:,i), colors(4));
end
%boxplot(Hyb_errors,'positions', agents,'symbol','', 'color',[0 0.5 0]);
%hold on
%boxplot(ICI_errors,'positions', agents,'symbol','','PlotStyle','compact')
%boxes = findall(gca,'Tag','Box');
%xticks(agents);
%xticklabels({'1','2','3','4','5','6','7','8','9','10'}) %,'80','100','120'
%hLegend = legend(boxes([1 end]), {'ICI','Hybrid'},'Fontsize',14);
legend('ground truth', 'Hybrid', 'ICI');
xlabel('Time Steps','FontSize', 12);
ylabel('State Values','Interpreter','latex','FontSize', 14);
%ax = gca;
%ax.FontSize = 12;
%ylim([0 1.1])


set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters 'cubesat_test_3.pdf'
%saveas(fig,'/home/naveed/Documents/DSE_data/Bhattacharyya_dist_10to100.pdf')
end