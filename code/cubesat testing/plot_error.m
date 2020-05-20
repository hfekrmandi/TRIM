function plot_error(error_)
global opt_dist
% 'Error' is Hyb/ICI. Row 1 is ICI, Row 2 is Hybrid
t = 1:opt_dist.i_step;
plot(getfield(error_(1), 'e_vs_cen'), t)