function nature_of_ci_gold_vs_silver(nAgent)
close all
clc
addpath(genpath(pwd))
eps_ = eps;
% addpath('C:\Users\Amirhossein\Dropbox\icra_2016_decentralized_est\code\problem_nature\graph_generators\randRegGraph')
% addpath('C:\Users\Amirhossein\Dropbox\icra_2016_decentralized_est\code\problem_nature\graph_generators\GrTheory')
% addpath('C:\Users\Amirhossein\Dropbox\icra_2016_decentralized_est\code\problem_nature\graph_generators\matlab_networks_routines\code')

% rand('state',0)
% test to see what happens if we change the topology (the network remains connected)
% nAgent = 10;
dim = 2;
global S1
% generate random covariance matrices
P_matrix = rand_cov_mat_gen_dim(dim*nAgent);
P_gold = inv(repmat(eye(dim),1,nAgent)*inv(P_matrix)*repmat(eye(dim),nAgent,1));

for i_agent = 1:nAgent
    range_ = (i_agent - 1)*dim + 1:(i_agent)*dim;
    S1(:,:,i_agent) = P_matrix(range_,range_);
end
P_GCI = global_ci(S1,'det',eps_);
figure
subplot(211)
for j_agent=1:nAgent
 h(j_agent)=error_ellipse(S1(:,:,j_agent));hold on;
 set(h(j_agent),'LineStyle','--')
end
 
% plot Gold
h(nAgent+1)=error_ellipse(P_gold);hold on;
 set(h(nAgent+1),'LineWidth',2)
% plot silver
h(nAgent+2)=error_ellipse(P_GCI);hold on;
 set(h(nAgent+2),'LineWidth',3)
subplot(212)
for j_agent=1:nAgent
 hh(j_agent)=error_ellipse(S1(:,:,j_agent));hold on;
 set(hh(j_agent),'LineStyle','--')
end
% legend('1','2','3','4','5')

% 
% % generate random topology
% for i=1:10
%     adj_tri = createRandRegGraph(nAgent, deg_con);
% %     bg2=biograph(adj_tri);
% %     bg2.layoutType = 'equilibrium';
% %     bg2.showarrows = 'off';
%     %     view(bg2);
%     bg2 =[];
%     G = generate_graph(full(adj_tri)+eye(nAgent));
%     if checkc(adj_tri)
%         result = iterative_ci(G,'det',eps_);
%         [h_fig,flag_disag] = plot_nature_of_ci(result,P_GCI,nAgent,adj_tri,bg2);
%         if flag_disag
%             disp('one_case_found')
% %             savefig(h_fig,['res1_nagent = ',num2str(nAgent),'_ deg_con = ',num2str(deg_con),'_ eps = ',num2str(eps_)]);
%             saveas(h,h_fig,['res1_nagent = ',num2str(nAgent),'_ deg_con = ',num2str(deg_con),'_ eps = ',num2str(eps_)]); 
%         end
%     end
%     
% end