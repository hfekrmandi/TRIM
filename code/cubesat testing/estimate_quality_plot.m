load('/home/naveed/Documents/DSE_data/Estimate_quality_10to100.mat') 
agents_array = [10 30 60 100];
N = size(agents_array,2);
range_step = 20;
index = 1;
 e_BC_dist_gold_Hybrid = zeros(range_step,N);
 e_BC_dist_gold_ICI = zeros(range_step,N);
 e_BC_dist_gold_gold = zeros(range_step,N);
 

            
for n_agents = agents_array
    %for i_agents = 1:n_agents
    for i_step = 1:range_step
        
        temp = x_gold{N,i_step,:};
        x_gold_100 = mean(temp);

        temp_mat = zeros(80,80);
        for i_agents = 1:100
            temp_mat = temp_mat + P_gold{N,i_step,i_agents};
        end
        temp_mat = temp_mat./100; 
        P_gold_100 = temp_mat;

        if n_agents == 100
            
            for i_agents = 1:n_agents
                
                temp_hyb(i_agents) = BC_distance(x_gold{N,i_step,i_agents},...
                        P_gold{N,i_step,i_agents},x_Hybrid{index,i_step,i_agents},P_Hybrid{index,i_step,i_agents});
                temp_ici(i_agents) = BC_distance(x_gold{N,i_step,i_agents},...
                        P_gold{N,i_step,i_agents},x_ICI{index,i_step,i_agents},P_ICI{index,i_step,i_agents});
                temp_gold(i_agents) = BC_distance(x_gold{N,i_step,i_agents},...
                        P_gold{N,i_step,i_agents},x_gold{index,i_step,i_agents},P_gold{index,i_step,i_agents});
            end
            
            e_BC_dist_gold_Hybrid(i_step,index) = mean(temp_hyb);
            e_BC_dist_gold_ICI(i_step,index) = mean(temp_ici);
            e_BC_dist_gold_gold(i_step,index) = mean(temp_gold);
        else
            

            temp = x_Hybrid{index,i_step,:};
            x_Hybrid_n_agents = mean(temp);

            temp = x_ICI{index,i_step,:};
            x_ICI_n_agents = mean(temp);

            temp = x_gold{index,i_step,:};
            x_gold_n_agents = mean(temp);

            temp_mat_hyb = zeros(80,80);
            temp_mat_ici = zeros(80,80);
            temp_mat_gold = zeros(80,80);

            for i_agents=1:n_agents
                temp_mat_hyb = temp_mat_hyb + P_Hybrid{index,i_step,i_agents};
                temp_mat_ici = temp_mat_ici + P_ICI{index,i_step,i_agents};
                temp_mat_gold = temp_mat_gold + P_gold{index,i_step,i_agents};
            end

            P_Hybrid_n_agents = temp_mat_hyb./n_agents;
            P_ICI_n_agents = temp_mat_ici./n_agents;
            P_gold_n_agents = temp_mat_gold./n_agents;

            temp = BC_distance(x_gold_100,P_gold_100, x_Hybrid_n_agents,P_Hybrid_n_agents);                                
            e_BC_dist_gold_Hybrid(i_step,index) = temp(1,1);

            temp = BC_distance(x_gold_100,P_gold_100, x_ICI_n_agents,P_ICI_n_agents);

            e_BC_dist_gold_ICI(i_step,index) = temp(1,1);

            temp = BC_distance(x_gold_100,P_gold_100, x_gold_n_agents,P_gold_n_agents);
            e_BC_dist_gold_gold(i_step,index) = temp(1,1);
            
            %trace
            trace_gold_100= trace(P_gold_100);
            trace_hybrid(i_step,index) = trace(P_Hybrid_n_agents)/trace_gold_100;
            trace_ICI(i_step,index) = trace(P_ICI_n_agents)/trace_gold_100;
            trace_gold(i_step,index) = trace(P_gold_n_agents)/trace_gold_100;
            
        end
    end
        
    %end
    index = index + 1;
end

%% plotting

fig = figure(1);
boxplot(e_BC_dist_gold_Hybrid, agents_array,'symbol','', 'color',[0 0.5 0])
hold on 
boxplot(e_BC_dist_gold_ICI, agents_array,'symbol','','PlotStyle','compact')
hold on
boxplot(e_BC_dist_gold_gold, agents_array,'symbol','', 'color',[0.9290, 0.6940, 0.1250])

xlabel('Number of receptors or agents')
ylabel('Bhattacharyya Distance')

boxes = findall(gca,'Tag','Box');
hLegend = legend(boxes([1 6 end]), {'FHS','ICI','Hybrid'},'Location','northwest');
xlabel('Number of receptors or agents')
ylabel('Bhattacharyya distance')
ylim([-0.09 1.1])

set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print -dpdf -painters '/home/naveed/Documents/DSE_data/estimate_quality_BC_distance_100.pdf'

%% 
%{
fig = figure(2);
plot([0,110],[1,1],'HandleVisibility','off')
hold on 

plot(agents_array,trace_ICI(range_step,:),'ob','MarkerFaceColor','b','MarkerSize',10)
hold on
plot(agents_array,trace_gold(range_step,:),'o','color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerSize',10)
hold on
plot(agents_array,trace_hybrid(range_step,:),'*','color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',10)
hold on 
xlim([0,110])
ylim([0.8,3])
legend('ICI','FHS','Hybrid')
xlabel('Number of receptors or agents')
ylabel('Average Trace of the Covariance Matrix (Normalised)')
set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters '/home/naveed/Documents/DSE_data/estimate_quality_trace.pdf'


%}























