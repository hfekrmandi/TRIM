function [x_next,P_next] = f(x_current,P_current,flag_noise)
global opt_dist
if flag_noise
    delta_noise =  randn(size(x_current))*sqrt(opt_dist.motion.Q);
    x_next = opt_dist.A*x_current + opt_dist.B*( opt_dist.source.Q') + delta_noise;
    P_next = opt_dist.A*P_current*opt_dist.A' + (opt_dist.motion.Q)*eye(size(x_next,1));
else
    x_next = opt_dist.A*x_current + opt_dist.B* opt_dist.source.Q';
    P_next = opt_dist.A*P_current*opt_dist.A'+(opt_dist.motion.Q)*eye(size(x_next,1));
end
end