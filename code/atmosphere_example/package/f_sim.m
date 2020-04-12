function x_next = f_sim(x_current,flag_noise)
global opt_dist
if flag_noise
    delta_noise =  randn(size(x_current))*sqrt(opt_dist.motion.Q);
    x_next = opt_dist.A*x_current + opt_dist.B*( opt_dist.source.Q') + delta_noise;
else
    x_next = opt_dist.A*x_current + opt_dist.B*opt_dist.source.Q';
end
end
