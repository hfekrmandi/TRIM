function covariance_intersection_test2()
clear all
clear global
global S1
dbstop if error
global nCovSamples;
nCovSamples= 50;
DegConnect = 4;
adj_mat = zeros(nCovSamples,nCovSamples);

% adj_mat = makerandCIJ_und(nCovSamples,ceil(nCovSamples*DegConnect));

%
%
for i=1:nCovSamples
    adj_mat(i,max(1,i-3):min(nCovSamples,i+3)) =1;
end

% adj_mat = [1 1 0 0;
%            1 1 1 0;
%            0 1 1 1;
%            0 0 1 1] ;

close all
% Sigma = [1.7 -1.5; -1.5 1.7];

df = 10;
for i=1:nCovSamples
    if i<=ceil(nCovSamples/2)
        Sigma = [1.7 -1.5; -1.5 1.7];
    else
        Sigma = [1 -.2; -.2 1];
        
    end
    S1(:,:,i) = wishrnd(Sigma,df)/df;
    h=error_ellipse(S1(:,:,i));
    set(h,'LineStyle','--');
    hold on
end
%% global CI
P_CI_tr = global_ci(S1,'tr');
P_CI_det = global_ci(S1,'det');

h_P_CI_tr=error_ellipse(P_CI_tr);

set(h_P_CI_tr,'Color',[1 0 0]);
set(h_P_CI_tr,'LineWidth',1.5);
hold on
h_P_CI_det=error_ellipse(P_CI_det);
set(h_P_CI_det,'Color',[0 1 0]);
set(h_P_CI_det,'LineWidth',1.5);
hold on

legend([h_P_CI_tr h_P_CI_det],'P (using trace)','P (using det)','Location','northoutside','Orientation','horizontal')


G = generate_graph(adj_mat);

%% iterative CI
result_ici_tr = iterative_ci(G,'tr');
result_ici_det = iterative_ci(G,'det');


%% MH
result_mh = mh_ci(G);


%% result comparison
figure
for i_cons=1:70
    i_cons
    h=error_ellipse(P_CI_tr);
    set(h,'Color',[1 0 0]);
    set(h,'LineWidth',1.5);
    
    h=error_ellipse(P_CI_det);
    set(h,'Color',[0 1 0]);
    set(h,'LineWidth',1.5);
    
    
    
    h=error_ellipse(inv(result_mh.consenus{1}.Y_prior{i_cons}));
    set(h,'Color',[0 1 0]);
    set(h,'LineWidth',1);
    hold on
    h=error_ellipse((result_ici_tr.consenus{1}.P_prior{i_cons}));
    set(h,'Color',[0 0 0]);
    set(h,'LineWidth',1);
    hold on
    h=error_ellipse((result_ici_det.consenus{1}.P_prior{i_cons}));
    set(h,'Color',[0 0 1]);
    set(h,'LineWidth',1);
    hold on
    % pause(0.1)
    hold off
end

plot_convergence_results(result_mh,'MH')
plot_convergence_results(result_ici_tr,'ICI_trace')
plot_convergence_results(result_ici_det,'ICI_det')


% % % % plot results
figure
subplot(121)
i_cons =70;
h_P_CI_tr=error_ellipse(P_CI_tr);
set(h_P_CI_tr,'Color',[1 0 0]);
set(h_P_CI_tr,'LineWidth',1.5);
hold on
h_P_CI_det=error_ellipse(P_CI_det);
set(h_P_CI_det,'Color',[0 1 0]);
set(h_P_CI_det,'LineWidth',1.5);

hold on
h_mh=error_ellipse(inv(result_mh.consenus{1}.Y_prior{i_cons}));
set(h_mh,'Color',[0 1 1]);
set(h_mh,'LineWidth',4);
hold on
h_ici_tr=error_ellipse((result_ici_tr.consenus{1}.P_prior{i_cons}));
set(h_ici_tr,'Color',[0 0 0]);
set(h_ici_tr,'LineWidth',1);
hold on

h_result_ici_det=error_ellipse((result_ici_det.consenus{1}.P_prior{i_cons}));
set(h_result_ici_det,'Color',[0 0 1]);
set(h_result_ici_det,'LineWidth',1);
hold on
for i=1:nCovSamples
    h=error_ellipse(S1(:,:,i));
    set(h,'LineStyle','--');
    hold on
end
axis tight
axis equal
legend([h_P_CI_tr h_P_CI_det,h_mh,h_ici_tr,h_result_ici_det],'global CI using trace',...
    'global CI using det','Local MH','Local ICI using trace','Local ICI using det')

subplot(122)
graph_draw(adj_mat)
axis tight
axis equal


end
% function cost_tr = cost_ci_tr(x)
% global S1_local
% information_matrix = zeros(size(S1_local(:,:,1)));
% for i=1:length(x)
%     information_matrix = information_matrix +x(i,1)*inv(S1_local(:,:,i));
% end
% cost_tr = trace(inv(information_matrix));
% end

% function cost_tr = cost_ci_det(x)
% global S1_local
% information_matrix = zeros(size(S1_local(:,:,1)));
% for i=1:length(x)
%     information_matrix = information_matrix +x(i,1)*inv(S1_local(:,:,i));
% end
% cost_tr = det(inv(information_matrix));
% end


% % function result = mh_ci(G)
% % global S1
% % 
% % for i_agent=1:size(G.Adj,2)
% %     result.consenus{i_agent}.Y_prior{1} = inv(S1(:,:,i_agent));
% % end
% % 
% % for i_consensus = 2:70
% %     for j_agent = 1 : size(G.Adj,2)
% %         result.consenus{j_agent}.Y_prior{i_consensus} = zeros(size(S1(:,:,1)));
% %         for k_agent = 1 :  size(G.Adj,2)
% %             p_jk = G.p(j_agent,k_agent);
% %             %             b_jk = covariance_intersection()
% %             updated_Y_prior = result.consenus{j_agent}.Y_prior{i_consensus} + p_jk*result.consenus{k_agent}.Y_prior{i_consensus-1};
% %             result.consenus{j_agent}.Y_prior{i_consensus} = updated_Y_prior;
% %         end
% %     end
% % end
% % end
