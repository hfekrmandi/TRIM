function nature_of_ci(nAgent,deg_con,eps_)
close all
clc
addpath(genpath(pwd))
% addpath('C:\Users\Amirhossein\Dropbox\icra_2016_decentralized_est\code\problem_nature\graph_generators\randRegGraph')
% addpath('C:\Users\Amirhossein\Dropbox\icra_2016_decentralized_est\code\problem_nature\graph_generators\GrTheory')
% addpath('C:\Users\Amirhossein\Dropbox\icra_2016_decentralized_est\code\problem_nature\graph_generators\matlab_networks_routines\code')

% rand('state',0)
% test to see what happens if we change the topology (the network remains connected)
% nAgent = 10;
dim = 2;
global S1
% generate random covariance matrices
for i_agent = 1:nAgent
    S1(:,:,i_agent) = rand_cov_mat_gen_dim(dim);
end
P_GCI = global_ci(S1,'det',eps_);

% generate random topology
for i=1:10
    adj_tri = createRandRegGraph(nAgent, deg_con);
%     bg2=biograph(adj_tri);
%     bg2.layoutType = 'equilibrium';
%     bg2.showarrows = 'off';
    %     view(bg2);
    bg2 =[];
    G = generate_graph(full(adj_tri)+eye(nAgent));
    if checkc(adj_tri)
        result = iterative_ci(G,'det',eps_);
        [h_fig,flag_disag] = plot_nature_of_ci(result,P_GCI,nAgent,adj_tri,bg2);
        if flag_disag
            disp('one_case_found')
%             savefig(h_fig,['res1_nagent = ',num2str(nAgent),'_ deg_con = ',num2str(deg_con),'_ eps = ',num2str(eps_)]);
            saveas(h,h_fig,['res1_nagent = ',num2str(nAgent),'_ deg_con = ',num2str(deg_con),'_ eps = ',num2str(eps_)]); 
        end
    end
    
end