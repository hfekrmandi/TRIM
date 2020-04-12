function sim_disconnect_gold()
global opt_dist fail_prob 
[flag_coon,G] = generate_graph_for_diag(opt_dist.reg_degree,opt_dist.n_degree_graph,fail_prob);
opt_dist.Graphs.G_obs = G;
opt_dist.dataG.Graph = G;
opt_dist.dataG.is_connected = flag_coon;
opt_dist.Graph_History{opt_dist.i_step} = opt_dist.dataG.Graph.p; 
end