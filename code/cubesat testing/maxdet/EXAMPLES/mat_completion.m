function A = mat_completion(Ainit,indices,C)
% A = mat_completion(Ainit,indices,C);
%
% Given A(x) = Ainit + \sum_{i=1}^m x_i (E_{r_i,c_i} + E_{c_i,r_i}) 
% and C=C^T, solve
%
%         minimize    Tr CA(x) - log det A(x)
%         subject to  A(x) > 0
%
% Ainit:    n-by-n positive definite matrix, gives the fixed part of A.
% indices:  2-by-m matrix, specifies the positions of the adjustable 
%           entries in the lower triangular part of A.
%           indices = [row entries; 
%                      col entries].
% C:        (optional) n-by-n symmetric matrix. If not given, C = 0.
%
% NOTE: mat_completion gives a demo if no input argument is given:

% last update: 03/19/96

% demo?
if nargin<1
  Ainit = eye(5)
  C = .1*ones(5,5)
  indices = [2 3 4 5 5 5 5;
             1 1 1 1 2 3 4]
elseif nargin<3
  C=zeros(size(Ainit));
end

% check input arguments
n = size(Ainit,1);
if (n ~= size(Ainit,2))
  error('Ainit must be square.');
end
Ainit=.5*(Ainit+Ainit');
if min(eig(Ainit)<=0)
  error('Ainit must be positive definite.')
end
if size(indices,1)~=2
  error('Indices must have two rows.');
end
row = indices(1,:);
col = indices(2,:);
m = length(row);
if (max(col-row)>=0)
  error('Indices must be in the strict lower triangular part of A.')
end

% form G(x)
G = [Ainit(:) zeros(n*n,m)];                % G_0
G((row-1)*n+col,[1:m]+1) = eye(m);
G((col-1)*n+row,[1:m]+1) = eye(m);          % G_i
G_blkszs = n;

% form F(x)  (F(x)=1 because the problem is unconstrained)
F = zeros(1,m+1);                           % F_i
F(1) = 1;                                   % F_0
F_blkszs = 1;

% form c (constant part of the objective is ignored in the optimization)
c = G(:,2:m+1)'*C(:);                       % c_i = Tr G_i*C

% call maxdet
x0 = zeros(m,1);                            % x0 = 0 is feasible
[x,Z,W,ul,hist,infostr] = maxdet(F,F_blkszs,G,G_blkszs,c,x0,...
      zeros(size(F,1),1),zeros(size(G,1),1),1e-8,0.5,100,100);

disp(infostr)

% form the solution A
A = reshape(G*[1;x],n,n);
