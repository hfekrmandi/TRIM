function [x,Z,W,ul,infostr] = phase1(F,F_blkszs,G,G_blkszs,gam,abstol,...
                                     reltol,NTiters);

% [x,Z,W,ul,infostr] = phase1(F,F_blkszs,G,G_blkszs,gam,abstol,reltol,...
%                             NTiters);
%
% Find an x s.t. F(x) > 0 and G(x) > 0
% or prove that no such x exists.
%
% phase1.m determines whether the optimal value of the semidefinite
% program
%  
%  minimize    t 
%  subject to  F(x) + t*I >=0 
%              G(x) + t*I >=0
% 
% is negative or not.
% 
% It produces either x and t with t<0.0 or a solution Z,W to 
% the dual problem
%
%  maximize    - Tr F_0*Z - Tr G_0*W
%  subject to  0 = Tr F_i*Z + Tr G_i*W, i=1,...,m
%              1 = Tr Z + Tr W
%              Z >= 0, W >= 0
% 
% with objective value non-negative.
%
% Input arguments:
% - F:        Block diagonal matrices F_i, i=0,...,m. 
%             F = [ F_0^1(:)  F_1^1(:) ...  F_m^1(:) ]
%                 [ F_0^2(:)  F_1^2(:) ...  F_m^2(:) ]
%                 [   ...       ...          ...     ]
%                 [ F_0^L(:)  F_1^L(:) ...  F_m^L(:) ]
%             F_i^j is the jth diagonal block of F_i.
%             It is assumed that the matrices diag(F_i,G_i) are linearly 
%             independent. 
%             WARNING: F will be modified if maxdet is terminated by an
%                      interrupt. 
% - F_blkszs: Integer L-vector. F_blkszs[j], j=0,...,L-1 gives the 
%             size of block j, ie, F_i^j has size F_blkszs[j] 
%             times F_blkszs[j].
% - G:        Block diagonal matrices G_i, i=0,...,m. 
%             G = [ G_0^1(:)  G_1^1(:) ...  G_m^1(:) ]
%                 [ G_0^2(:)  G_1^2(:) ...  G_m^2(:) ]
%                 [   ...       ...          ...     ]
%                 [ G_0^K(:)  G_1^K(:) ...  G_m^K(:) ]
%             G_i^j is the jth diagonal block of G_i.
%             It is assumed that the matrices diag(F_i,G_i) are linearly 
%             independent. 
%             WARNING: G will be modified if maxdet is terminated by an
%                      interrupt. 
% - G_blkszs: Integer K-vector. G_blkszs[j], j=0,...,K-1 gives the
%             size of block j, ie, G_i^j has size G_blkszs[j]
% - gam:      > 0. Controls how aggressively the algorithm works.
% - abstol:   absolute tolerance, >= MINABSTOL (MINABSTOL given in maxdet.h).
% - reltol:   relative tolerance, >=0.
% - NTiters:  maximum number of Newton (and predictor) iterations, >=1
%
% Output arguments:
% x:         m-vector; last primal iterate
% Z:         last dual iterate;  block-diagonal matrix stored as 
%            [ Z^1(:); Z^2(:); ... ; Z^L(:) ]
%            where Z^i is the ith diagonal block.
% W:         last dual iterate;  block-diagonal matrix stored as 
%            [ W^1(:); W^2(:); ... ; W^L(:) ]
%            where Z^i is the ith diagonal block.
% ul:        2-vector, primal and dual objective.
%            If ul(1) < 0.0, the problem is feasible and x is a 
%            solution.  If ul(2) > 0.0, the problem is infeasible
%            and Z,W form a proof of infeasibility.  If ul(1) > 0.0
%            and ul(2) < 0.0, no conclusion can be made.
% infostr:   'feasible' or
%            'infeasible' or
%            'feasibility cannot be determined'

m = size(F,2)-1;
if (m ~= size(G,2)-1)
  error('F and G must have the same number of columns.');
end
if (size(F,1) ~= sum(F_blkszs.^2)) 
  error('Dimensions of F do not match F_blkszs.');
end;
if (size(G,1) ~= sum(G_blkszs.^2)) 
  error('Dimensions of G do not match G_blkszs.');
end;

% mineigF is the smallest eigenvalue of F_0
mineigF = 0.0;
k=0; for n=F_blkszs,
   mineigF = min(mineigF, min(eig(reshape(F(k+[1:n*n],1),n,n))));
   k=k+n*n;   % k = sum n_i*n_i 
end;
% mineigG is the smallest eigenvalue of G_0
mineigG = 0.0;
k=0; for n=G_blkszs,
   mineigG = min(mineigG, min(eig(reshape(G(k+[1:n*n],1),n,n))));
   k=k+n*n;   % k = sum n_i*n_i 
end;

% eyeF is the identity
eyeF = zeros(size(F,1),1);  
k=0; for n=F_blkszs,
   eyeF(k+[1:n*n]) = reshape(eye(n),n*n,1);   % identity
   k=k+n*n;   % k = sum n_i*n_i 
end;
% eyeG is the identity
eyeG = zeros(size(G,1),1);  
k=0; for n=G_blkszs,
   eyeG(k+[1:n*n]) = reshape(eye(n),n*n,1);   % identity
   k=k+n*n;   % k = sum n_i*n_i 
end;

% initial x0 
x0 = [zeros(m,1); max(-1.1*min(mineigF,mineigG), 1)];

% linear objective
c = [zeros(m,1); 1];

% call maxdet
[x,Z,W,ul,hist,infostr]=maxdet([F,eyeF; G,eyeG],...
                               [F_blkszs(:)',G_blkszs(:)'],...
                               [1 zeros(1,m+1)],1,c,x0,...
                               zeros(size(F,1)+size(G,1),1),0,...
                               abstol,reltol,gam,NTiters);

% prepare output
x = x(1:m);
W = Z(size(F,1)+1:size(F,1)+size(G,1));
Z = Z(1:size(F,1));
if ul(1)<0
  infostr = 'feasible';
elseif ul(2)>=0
  infostr = 'infeasible';
else
  infostr = 'feasibility cannot be determined';
end