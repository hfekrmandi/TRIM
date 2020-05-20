clearvars
close all
clc
load('80_states_convergence_rate_40_60.mat')
agents_count = n_receptors_array;
steps = size(converg_steps_array,2);
T = size(error_(1,1).e_BC_dist_gold,2);    
N = size(agents_count,2);

Hyb_errors = zeros(T, N,steps);
ICI_errors = zeros(T, N,steps);
time = zeros(N,steps);
fig = figure(1);

for i=[2 1]
    for j=1:steps
        Hyb_errors(:,i,j) = error_(i,j).e_BC_dist_gold(1,:)./mean(error_(i,end).e_BC_dist_gold(1,:));
        %ICI_errors(:,i,j) = error_(i,j).e_BC_dist_gold(2,:)./mean(error_(i,end).e_BC_dist_gold(2,:));
        %time(i,j) = time_array(i,j);
    end


    
    xlim = [0 210];
    if i == 1
        color = [1 0 0];
    else
        color = [0 0 1];
    end
    boxplot(squeeze(Hyb_errors(:,i,:)),'positions', converg_steps_array,'symbol','', 'color',color)
    
    hold on 
    %boxplot(squeeze(ICI_errors(:,i,:)),'positions', converg_steps_array,'symbol','','PlotStyle','compact')
   
    
   
    
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
hLegend = legend(boxes([1 end]), {'Hybrid 40 agents','Hybrid 60 agents'},'Location','northwest','Fontsize',14);
xlabel('Number of consensus iterations')
ylabel('Bhattacharyya distance (normalised)')
ylim([0 1.2])

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print -dpdf -painters 'convergence_80states_40_60.pdf' ;
%print -dpdf -painters 'time_taken_80states_30agents.pdf' ;