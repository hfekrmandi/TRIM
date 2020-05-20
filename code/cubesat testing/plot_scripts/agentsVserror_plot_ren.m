clearvars
close all
clc
load('/home/naveed/Documents/DSE_data/80_states_10to100_Ren.mat')
agents_count = n_receptors_array;

T = size(error_(1).e_BC_dist_gold,2);    
N = size(agents_count,2);
Hyb_errors_ren = zeros(T, N);
ICI_errors = zeros(T, N);
for i=1:N
    Hyb_errors_ren(:,i) = error_(i).e_BC_dist_gold(1,:);
    ICI_errors(:,i) = error_(i).e_BC_dist_gold(2,:);
end

load('/home/naveed/Documents/DSE_data/80_states_10to100.mat')
agents_count = n_receptors_array;

T = size(error_(1).e_BC_dist_gold,2);    
N = size(agents_count,2);
Hyb_errors = zeros(T, N);
%ICI_errors = zeros(T, N);
for i=1:N
    Hyb_errors(:,i) = error_(i).e_BC_dist_gold(1,:);
    %ICI_errors(:,i) = error_(i).e_BC_dist_gold(2,:);
end

fig = figure(1);
boxplot(Hyb_errors_ren,'positions', agents_count,'symbol','', 'color',[1 0 0])
hold on 
boxplot(ICI_errors,'positions', agents_count,'symbol','','PlotStyle','compact')
hold on 
boxplot(Hyb_errors,'positions', agents_count,'symbol','', 'color',[0 0.5 0])
hold on
boxes = findall(gca,'Tag','Box');
xticks(agents_count)
xticklabels({'10','15','20','25','30','40','50','60','70','80','90','100'}) 
hLegend = legend(boxes([1 end 13]), {'Hybrid','Hybrid 1-step','ICI'},'Fontsize',14);
xlabel('Number of sensors','FontSize', 12)
ylabel('$D_B$','Interpreter','latex','FontSize', 14)
ax = gca;
ax.FontSize = 12; 

ylim([0 1.1])

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters '/home/naveed/Dropbox/Research/Data/T_RO_DSE/Bhattacharyya_dist_80state_10to100_Ren.pdf'
%saveas(fig,'/home/naveed/Documents/DSE_data/Bhattacharyya_dist_10to100.pdf')