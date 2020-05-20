function [lambda,S] = exp_design(V)
% [lambda,S] = exp_design(V);
%
% D-optimal experiment design:
%   minimize   -log det sum_i lambda_i v_iv_i^T
%   subject to  lambda_i >= 0
%               sum_i lambda_i == 1
%
% V:      pxM matrix, V=[v_1 v_2 ... v_M], contains M test vectors 
%         in R^p
% lambda: lambda_i is the fraction of the experiments allocated to
%         test vector v_i
% S:      sum_i lambda_i v_iv_i^T
%

% last update: 02/20/96

[n,m]=size(V);
mm=m-1;    % number of variables: [\lambda_1 ... \lambda_{m-1}]

% form G
G=zeros(n*n,mm+1);
tmpmat=V(:,m)*V(:,m)';
G(:,1)=reshape(tmpmat,n*n,1);
for i=1:m-1
  G(:,i+1)=reshape(V(:,i)*V(:,i)'-tmpmat,n*n,1);
end
G_blkszs=n;

% form F
F=[];  F_blkszs=[];
% non-negativity on lambda_i
for i=1:m-1
  F=[F; [0 zeros(1,i-1) 1 zeros(1,m-1-i)]];
  F_blkszs=[F_blkszs 1];
end
F=[F; [1 -ones(1,m-1)]];
F_blkszs=[F_blkszs 1];

% lambda_i = 1/m is initial feasible
x0 = (1/m)*ones(mm,1);    

% call maxdet
c = zeros(mm,1);
[x,Z,W,ul,hist,infostr]=maxdet(F,F_blkszs,G,G_blkszs,c,x0,...
                               zeros(size(F,1),1),zeros(size(G,1),1),...
                               1e-5,1e-5,100,100);

% form lambda
lambda=x(1:m-1);
lambda=[lambda; 1-sum(lambda)];

% form (1/N)A^TA
S=lambda(1)*V(:,1)*V(:,1)';
for i=2:m
  S=S+lambda(i)*V(:,i)*V(:,i)';
end
