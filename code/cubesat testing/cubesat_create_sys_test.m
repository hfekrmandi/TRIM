function [A,x0,B,C] = cubesat_create_sys_test(n, p, m, nAgents)
global opt_dist

% Creates a block diagonal matrix of n-th order models each with p outputs and m inputs
% This model keeps each bot stationary, only effected by the input
A0 = eye(n);
B0 = eye([n, m]);
C0 = eye([p, n]);
D0 = zeros([p, m]);

for i = 1:nAgents
    A_list{i} = A0;
    B_list{i} = B0;
    C_list{i} = C0;
    D_list{i} = D0;
    source.Q(i*n - n + 1:i*n) = i:i+n-1;
end

A = blkdiag(A_list{:});
B = blkdiag(B_list{:});
C = blkdiag(C_list{:});
D = blkdiag(D_list{:});

sys = ss(A, B, C, D);
rank(obsv(sys.A,sys.C))
rank(ctrb(sys.A,sys.B))

opt_dist.sys = sys;
opt_dist.A = full(sys.A);
opt_dist.B = full(sys.B);
opt_dist.C = full(sys.C);

x0 = zeros(size(A,1),1);
opt_dist.source = source;

opt_dist.nAgents = size(opt_dist.C,1);