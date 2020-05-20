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

for i=1
    for j=1:steps
        Hyb_errors(:,i,j) = error_(i,j).e_BC_dist_gold(1,:);
        ICI_errors(:,i,j) = error_(i,j).e_BC_dist_gold(2,:);
        %time(i,j) = time_array(i,j);
    end



    
    fig = figure(i);
    boxplot(squeeze(Hyb_errors(:,i,:)), converg_steps_array,'symbol','', 'color',[0 0.5 0])
    %hLegend = legend(findall(gca,'Tag','Box'), {'Hybrid'});
    hold on 
    boxplot(squeeze(ICI_errors(:,i,:)), converg_steps_array,'symbol','','PlotStyle','compact')
    boxes = findall(gca,'Tag','Box');
    hLegend = legend(boxes([1 end]), {'ICI','Hybrid'});
    xlabel('Number of consensus iterations')
    ylabel('Bhattacharyya distance')
    ylim([0 1.1])

    set(fig,'Units','inches');
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    
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
print -dpdf -painters '/home/naveed/Documents/DSE_data/Converg_rate_80states_30agents.pdf' ;





