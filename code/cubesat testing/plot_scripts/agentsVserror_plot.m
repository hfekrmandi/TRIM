function agentsVserror_plot()
clearvars
close all
clc
%{
load('80_states_10to100_Ren.mat')
agents_count = n_receptors_array;

T = size(error_(1).e_BC_dist_gold,2);    
N = size(agents_count,2);
Hyb_errors_ren = zeros(T, N);
ICI_errors = zeros(T, N);
for i=1:N
    Hyb_errors_ren(:,i) = error_(i).e_BC_dist_gold(1,:);
    %ICI_errors(:,i) = error_(i).e_BC_dist_gold(2,:);
end
%}
load('cubesat_test_3.mat')
N = opt_dist.nAgents;
agents = 1:N;

T = size(error_results{1, 1}.error_Hybrid.e_BC_dist_gold,2);
Hyb_errors = zeros(T, N);
ICI_errors = zeros(T, N);
for i=1:N
    Hyb_errors(:,i) = error_results{1, i}.error_Hybrid.e_BC_dist_gold(1,:);
    ICI_errors(:,i) = error_results{1, i}.error_ICI.e_BC_dist_gold(1,:);
end
%%
fig = figure(1);
boxplot(Hyb_errors,'positions', agents,'symbol','', 'color',[0 0.5 0])
%hLegend = legend(findall(gca,'Tag','Box'), {'Hybrid'});
hold on 
boxplot(ICI_errors,'positions', agents,'symbol','','PlotStyle','compact')
boxes = findall(gca,'Tag','Box');
xticks(agents)
xticklabels({'1','2','3','4','5','6','7','8','9','10'}) %,'80','100','120'
hLegend = legend(boxes([1 end]), {'ICI','Hybrid'},'Fontsize',14);
xlabel('Time Steps','FontSize', 12)
ylabel('$D_B$','Interpreter','latex','FontSize', 14)
ax = gca;
ax.FontSize = 12; 
%xl = xlabel(printnombrejpg, 'FontSize', 14);
%get(xl);
ylim([0 1.1])

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters 'cubesat_test_3.pdf'
%saveas(fig,'/home/naveed/Documents/DSE_data/Bhattacharyya_dist_10to100.pdf')
end