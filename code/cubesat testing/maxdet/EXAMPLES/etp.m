function [d,Q] = etp(A)
% [d,Q] = etp(A);
%
% solves the eductional testing problem:
% (primal problem) maximize     e'*d 
%                  subject to   A - diag(d) >= 0
%                               d >= 0
%
% (dual problem)   minimize     Tr A*Q
%                  subject to   Q >= 0
%                               diag(Q) >= e
% 
% (e is the vector with all components one)
% 
% - A: nxn, positive definite matrix.
% - d: n-vector, solution for primal problem.
% - Q: nxn, solution for dual problem.
%

% last update: 03/18/96

n = size(A,1);  
if (size(A,2) ~= n) error('A must be square.'); end;
if (norm(A - A') > 1e-10*norm(A)) error('A must be symmetric.'); end;
lambdamin = min(eig(A)); 
if lambdamin < 1e-8, error('A must be positive definite.');  end;

% c = -e
c = -ones(n,1);    

% F = [A -E1 ... -En]  (Ei: zero except for E_ii = 1) 
F = [A zeros(n,n*n)];
F(:,n+[1:n+1:n*n]) = -eye(n);
   
% F = [ A(:) -E1(:) ... -En(:); 0 e1 ... en ] (ei: i-th unit vector)
F = [reshape(F,n*n,n+1); zeros(n,1) eye(n)];
F_blkszs = [n; ones(n,1)];

% G = 1
G = zeros(1,size(F,2));
G(1) = 1;

% primal initial point
x0 = 0.9*lambdamin*ones(n,1); 

% call maxdet
[d,Z,W,ul,hist,infostr] = maxdet(F,F_blkszs,G,1,c,x0, ...
   zeros(size(F,1),1),zeros(size(G,1),1),1e-3,1e-3,100,100);
Q = reshape(Z(1:n*n),n,n);
