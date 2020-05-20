function generate_diagram_rel_perf()
clear all
close all
dbclear if error
% dbclear if warning

% figure

n_degree_graph = 15;
nreg =10;



% adj_tri = ones(n_degree_graph,n_degree_graph);
% graph_draw(adj_tri)
res_diag = 100;
nStep_test = 1000;
rel = zeros(res_diag,nStep_test);
idx_reg = 2:2:n_degree_graph;


for i_reg=1:length(idx_reg)
    eln = kregular(n_degree_graph,idx_reg(i_reg));
    adj_tri = zeros(n_degree_graph,n_degree_graph);
    idx = sub2ind(size(adj_tri), eln(:,1)', eln(:,2)');
    adj_tri(idx)=1;
    for i_res=1:res_diag
        thresh = i_res/res_diag;
        for i_step =1:nStep_test
            rand_mask = rand(n_degree_graph,n_degree_graph);
            
            adj_ = double(or(adj_tri.*rand_mask>=thresh,eye(n_degree_graph)));
            tri_u1 = ([triu(adj_,1)' + triu(adj_)]);
            adj_ = tri_u1;
            if(isconnected(adj_))
                rel(i_res,i_step,i_reg) = 1;
            else
                %         graph_draw(adj_)
            end
        end
    end
end
thresh_ = 1/res_diag:1/res_diag:1;
figure
for i_plot=1:size(rel,3)
    h(i_plot) = plot(thresh_,100*sum(rel(:,:,i_plot)/nStep_test,2),'DisplayName',['(n = ',num2str(n_degree_graph) ,', k = ',num2str(idx_reg(i_plot)),')' ],'LineWidth',2); hold on
    legend('-DynamicLegend');

%     legend( h(i_plot),['(n = ',num2str(n_degree_graph) ,', k = ',num2str(idx_reg(i_plot)),')' ])
end
xlabel('Probablity of Failure')
ylabel('Percentage of time Graph is connected')

