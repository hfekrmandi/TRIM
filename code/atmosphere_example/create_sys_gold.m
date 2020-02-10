function [A,x0,B,C] = create_sys_gold(n,p,m)
global opt_dist
% drss(n,p,m) generates an n-th order model with p outputs and m inputs.
is_stable = 0;
while ~is_stable
    sys = drss(n,p,m);
    is_stable = isstable(sys);
end
rank(obsv(sys.A,sys.C));
rank(ctrb(sys.A,sys.B));

A = full(sys.A);
B = full(sys.B);
C = full(sys.C);
opt_dist.sys = sys;
opt_dist.A = full(sys.A);
opt_dist.B = full(sys.B);
opt_dist.C = full(sys.C);

x0 =zeros(size(A,1),1);
source.Q = ones(1,m);
opt_dist.source = source;


