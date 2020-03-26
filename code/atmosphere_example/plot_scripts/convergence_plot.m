clearvars
close all
clc
load('/home/naveed/Documents/DSE_data/80_states_convergence_rate_30_100.mat')
agents_count = n_receptors_array;
steps = size(converg_steps_array,2);
T = size(error_(1,1).e_BC_dist_gold,2);    
N = size(agents_count,2);

Hyb_errors = zeros(T, N,steps);
ICI_errors = zeros(T, N,steps);
time = zeros(N,steps);
fig = figure(1);
for i=1 
    for j=1:steps
        Hyb_errors(:,i,j) = error_(i,j).e_BC_dist_gold(1,:);
        ICI_errors(:,i,j) = error_(i,j).e_BC_dist_gold(2,:);
        %time(i,j) = time_array(i,j);
    end


    
    xlim = [0 210];

    h = boxplot(squeeze(Hyb_errors(:,i,:)),'positions', converg_steps_array,'symbol','', 'color',[0 0.5 0]);
    set(h,'LineWidth',2)
    %hLegend = legend(findall(gca,'Tag','Box'), {'Hybrid'});
    hold on 
    boxplot(squeeze(ICI_errors(:,i,:)),'positions', converg_steps_array,'symbol','','PlotStyle','compact')
   
    
   
    
    %{
    fig = figure(i);
    plot(converg_steps_array,time(i,:),'*','color',[0 0.5 0])
    xlabel('Number of consensus iterations')
    ylabel('Time taken')
    set(fig,'Units','inches');
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    %}
   
end
%%

boxes = findall(gca,'Tag','Box');
xticks(converg_steps_array)
xticklabels({'1','','10','','20','30','40','50','60','80','100','120','140','160','180','200'}) %,
hLegend = legend(boxes([1 end]), {'ICI','Hybrid'},'Location','northwest','Fontsize',14);
xlabel('Number of consensus iterations','FontSize', 12)
ylabel('$D_B$','Interpreter','latex','FontSize', 14)
ax = gca;
ax.FontSize = 12; 
ylim([0 1.1])

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print -dpdf -painters '/home/naveed/Dropbox/Research/Data/T_RO_DSE/convergence_80states_30agents.pdf' ;
%print -dpdf -painters 'time_taken_80states_30agents.pdf' ;