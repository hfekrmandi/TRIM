clearvars
close all
clc
load('80_states_convergence_rate_40_60.mat')
agents_count = n_receptors_array;
steps = size(converg_steps_array,2);
T = size(error_(1,1).e_BC_dist_gold,2);    
N = size(agents_count,2);

Hyb_errors = zeros(N,steps);
ICI_errors = zeros(N,steps);
time = zeros(N,steps);
fig = figure(1);
%%
for i=[1 2]
    for j=1:steps
        Hyb_errors(i,j) = mean(error_(i,j).e_BC_dist_gold(1,:))/mean(error_(i,end).e_BC_dist_gold(1,:));
        %ICI_errors(:,i,j) = error_(i,j).e_BC_dist_gold(2,:)./mean(error_(i,end).e_BC_dist_gold(2,:));
        %time(i,j) = time_array(i,j);
    end


    
    xlim = [0 210];
    if i == 1
        color = [1 0 0];
    else
        color = [0 0 1];
    end
    
    plot(converg_steps_array, Hyb_errors(i,:), '-o','MarkerEdgeColor',color, 'MarkerSize',10, 'color',color); 
    hold on 
end
%%
load('80_states_convergence_rate_30_100.mat')
agents_count = n_receptors_array;
steps = size(converg_steps_array,2);
T = size(error_(1,1).e_BC_dist_gold,2);    
N = size(agents_count,2);

Hyb_errors = zeros(N,steps);
ICI_errors = zeros(N,steps);
time = zeros(N,steps);

for i=1
    for j=1:steps
        Hyb_errors(i,j) = mean(error_(i,j).e_BC_dist_gold(1,:))/mean(error_(i,end).e_BC_dist_gold(1,:));
        %ICI_errors(:,i,j) = error_(i,j).e_BC_dist_gold(2,:)./mean(error_(i,end).e_BC_dist_gold(2,:));
        %time(i,j) = time_array(i,j);
    end

    if i == 1
        color = [0 1 0];
    else
        color = [0 0 0];
    end
    
    plot(converg_steps_array, Hyb_errors(i,:), '-o','MarkerEdgeColor',color, 'MarkerSize',10, 'color',color); 
    hold on 
end

%%
converg_steps_array_old = converg_steps_array;
load('80_states_convergence_rate_10_40_100.mat')
agents_count = n_receptors_array;
steps_new = size(converg_steps_array,2);
%T = size(error_(1,1).e_BC_dist_gold,2);    
N = size(agents_count,2);

Hyb_errors = zeros(N,steps);
ICI_errors = zeros(N,steps);
time = zeros(N,steps);

for i=1
    for j=1:steps
        if j<=steps_new
            Hyb_errors(i,j) = mean(error_(i,j).e_BC_dist_gold(1,:))/mean(error_(i,end).e_BC_dist_gold(1,:));
        else
            Hyb_errors(i,j) = mean(error_(i,steps_new).e_BC_dist_gold(1,:))/mean(error_(i,end).e_BC_dist_gold(1,:));
        end  
    end
   
    plot(converg_steps_array_old, Hyb_errors(i,:), '-o','MarkerEdgeColor','m', 'MarkerSize',10, 'color','m'); 
    hold on 
end

%%
%boxes = findall(gca,'Tag','Box');
%xticks(converg_steps_array)
%xticklabels({'1','','10','','20','30','40','50','60','80','100','120','140','160','180','200'}) %,
hLegend = legend({'40 agents','60 agents','30 agents','10 agents'},'Location','southeast','Fontsize',14);

labels = get(legend(), 'String');
plots = flipud(get(gca, 'children'));

% Now re-create the legend
%neworder = [5, 3, 1, 2,4];
neworder = [4,3,1,2];
legend(plots(neworder), labels(neworder))

%neworder = [5, 3, 1, 2,4];
%hLegend.PlotChildren = hLegend.PlotChildren(neworder);

xlabel('Number of consensus iterations')

ylabel('Mean $D_B$ (normalised)','Interpreter','latex')
ylim([0 1.2])

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print -dpdf -painters 'convergence_80states_means_1.pdf' ;
%print -dpdf -painters 'time_taken_80states_30agents.pdf' ;