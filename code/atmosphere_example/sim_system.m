function sim_system()
global opt_dist
sim_disconnect()
opt_dist.sim.gt.x_bar = f_sim(opt_dist.result.gt.x_bar,1);
opt_dist.result.gt.x_bar = opt_dist.sim.gt.x_bar ;
sim_obs();