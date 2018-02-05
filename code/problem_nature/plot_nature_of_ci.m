function [h_fig,flag_disag]= plot_nature_of_ci(result,P_GCI,nAgent,adj_tri,bg2)
   Dia_graph =  diameter(adj_tri);
% persistent error_results
max_it = 70;
global S1
flag_disag = 0;
h_fig = -1;
for i_consensus=1:150;
    est{i_consensus}.P_cen = P_GCI;
       
    
    
    for j_agent = 1 :nAgent
        
        est{i_consensus}.P_bar(:,:,j_agent) = result.consenus{j_agent}.P_prior{i_consensus};
        e_P_bar_cen(i_consensus,j_agent) = norm(est{i_consensus}.P_bar(:,:,j_agent) - est{i_consensus}.P_cen );

    end
end
for j_agent = 1 :nAgent
    
    if  e_P_bar_cen(i_consensus,j_agent)>=0.01
        flag_disag = 1;
    end
end


error_results.e_P_bar_cen = e_P_bar_cen;




close all
if flag_disag
h_fig= figure;
subplot(311)
for j_agent=1:nAgent
plot( error_results.e_P_bar_cen(:,j_agent) ) ; hold on;
end
title('e_{P_{bar\_cen}} = norm(P_{bar} - P_{cen})')
legend('1','2','3','4','5')
subplot(312)
for j_agent=1:nAgent
 h(j_agent)=error_ellipse(result.consenus{j_agent}.P_prior{max_it});hold on;
end
 h(j_agent+1)=error_ellipse(P_GCI);
set(h(j_agent+1),'LineWidth',3)
for j_agent=1:nAgent
 h(nAgent+1+j_agent)=error_ellipse(S1(:,:,j_agent));hold on;
 set(h(nAgent+1+j_agent),'LineStyle','--')

end
title(['Diameter of the Graph = ',num2str(Dia_graph)]);
subplot(313)
graph_draw(adj_tri)
axis equal
end
