function problem_def_gold(A,B,C,x0,converg_steps)
global opt_dist

% log file 
opt_dist.fid_log = fopen('log.txt','at');
fprintf(opt_dist.fid_log , 'text line number 1 \n');

% opt_dist.save_dir = uigetdir();
% if set to 1 a set of sanity checks are performed duting the system setup.
% For example, rows of the G should add up to 1, the system should be 
% observable.
opt_dist.FLAGS.sanity_checks = 1; 
opt_dist.FLAGS.compare_with_CI = 1; % if set to 1
opt_dist.FLAGS.our_method = 1; % if set to 1
opt_dist.FLAGS.pure_CI = 0; % if set to 1
opt_dist.FLAGS.debug_CI = 0; % if set to 1
opt_dist.FLAGS.verbose = 0; % if set to 1
opt_dist.FLAGS.debug = 0;
opt_dist.FLAGS.debug_consensus = 0;
opt_dist.FLAGS.obs_noise_type =  'absolute'; % ('absolute' | 'relative');

opt_dist.figures.fig_cov_debug = figure;

opt_dist.nAgents = size(opt_dist.C,1); %Number of agents
opt_dist.dt = 1;
opt_dist.dimAgents = 3;
opt_dist.obs.Range = 1.2;
opt_dist.dimState = size(A,1); % 
q_var = 10^(-6);
r_var = 5*10^(-5);

% Q = q_var*eye(size(A));
% R = r_var*eye(size(C,1));
% For definition of regular graphs and 
% regularity degree look at https://en.wikipedia.org/wiki/Regular_graph
opt_dist.reg_degree = 4; 
opt_dist.n_degree_graph = opt_dist.nAgents;
opt_dist.obs.rel_perc = 0.05;
opt_dist.obs.R = r_var;
opt_dist.iter_interest = [1];
opt_dist.dimObs =1;%opt_dist.dimAgents;
opt_dist.nIterations = 60;
opt_dist.nSteps = converg_steps + 1; %consensus steps (change to 2 for Ren implementation) +1 for initialisation step
opt_dist.scenario = '2';
opt_dist.motion.Q = q_var;%(opt_dist.source.Q.*0.1).^2;



%%%% building the graph
full_Adj = ones(opt_dist.nAgents,opt_dist.nAgents);
topol_Adj = zeros(opt_dist.nAgents,opt_dist.nAgents);
for i=1:opt_dist.nAgents
    topol_Adj(i,max(1,i-3):min(opt_dist.nAgents,i+3)) =1;
end
obs_Adj = topol_Adj;
opt_dist.A = A;
G_full = generate_graph(full_Adj);
G = generate_graph(topol_Adj);
G_obs = generate_graph(obs_Adj);

opt_dist.Graphs.G_obs = G_obs;
opt_dist.Graphs.G_full = G_full;
opt_dist.Graphs.G = G;
nv = opt_dist.nAgents;
x_state(:) = zeros(1,opt_dist.dimState);
opt_dist.x_gt = x_state(:);
opt_dist.result.gt.x_bar = opt_dist.x_gt;

x0 = opt_dist.x_gt;

% figure
% imagesc(full(G_fault.Adj))
% imagesc(double(G.Adj))

if(opt_dist.FLAGS.sanity_checks)
    for i=1:nv
        disp([i,sum(G.p(i,:))]); %??
    end
end

i_time = 1;
opt_dist.i_time = i_time;
opt_dist.result.prior.x_cen =  x_state(:);
opt_dist.result.prior.P_cen =  0.05.*eye(opt_dist.dimState);
for i_agent = 1 : opt_dist.nAgents
    opt_dist.result.prior.x_bar(:,i_agent) = opt_dist.x_gt;% + randn(size(opt_dist.x_gt,1),1).*sqrt(0.05);
    opt_dist.result.prior.P_bar(:,:,i_agent) = 0.05*eye(opt_dist.dimState);
    if opt_dist.FLAGS.compare_with_CI
        opt_dist.result.prior.x_bar_CI(:,i_agent) = opt_dist.result.prior.x_bar(:,i_agent);
        opt_dist.result.prior.P_bar_CI(:,:,i_agent) = opt_dist.result.prior.P_bar(:,:,i_agent);
    end
    opt_dist.result.initial.x_bar(:,i_agent) = opt_dist.x_gt;% + randn(size(opt_dist.x_gt,1),1).*sqrt(0.05);
    opt_dist.result.initial.P_bar(:,:,i_agent) = 0.05*eye(opt_dist.dimState);
    if opt_dist.FLAGS.compare_with_CI
        opt_dist.result.initial.x_bar_CI(:,i_agent) = opt_dist.result.prior.x_bar(:,i_agent);
        opt_dist.result.initial.P_bar_CI(:,:,i_agent) = opt_dist.result.prior.P_bar(:,:,i_agent);
    end
end
