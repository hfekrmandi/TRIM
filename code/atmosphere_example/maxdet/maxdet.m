% [x,Z,W,ul,hist,infostr]=maxdet(F,F_blkszs,G,G_blkszs,c,x0,Z0,W0,abstol,...
%                                reltol,gam,NTiters);
%
% Solves determinant-maximization (MAXDET) program 
%
% (P) minimize    c^Tx - \log\det G(x)
%     subject to  G(x)\geq 0
%                 F(x)\geq 0
%
% and its dual
%
% (D) maximize    \log\det W -\Tr ZF_0 -\Tr WG_0 + l
%     subject to  W\geq 0, Z\geq 0
%                 \Tr ZF_i + \Tr WG_i = c_i, i=1,...,m
%
% Convergence criteria: 
% (1) maximum number of Newton iterations is exceeded
% (2) duality gap is less than abstol
% (3) primal and dual objective are both positive and 
%     duality gap is less than reltol * dual objective;         
%     or, primal and dual objective are both negative and
%     duality gap is less than reltol * minus the primal objective
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
% - c:        m-vector, specifies primal objective.
% - x0:       m-vector, strictly primal feasible point. 
% - Z0:       [ Z0^1(:); ... ; Z0^L(:) ], not used if NOT strictly dual
%             feasible.
% - W0:       [ W0^1(:); ... ; W0^K(:) ], not used if NOT strictly dual
%             feasible.
% - gam:      algorithm parameter which affects convergence rate, >0.
% - abstol:   absolute tolerance, >= MINABSTOL (MINABSTOL given in maxdet.h).
% - reltol:   relative tolerance, >= 0.
% - NTiters:  maximum total number of Newton (and predictor) iterations, >=1.
%
% Output arguments:
% - x, Z, W:  last primal and dual iterate.
% - ul:       two-vector.  ul(1) = c^Tx -\logdet G(x),
%             ul(2) = \logdet W -\Tr ZF_0 -\Tr WG_0 + l.
% - hist:     history. history=[objective;
%                               duality gap;
%                               number of Newton iterations]
% - infostr:  'maximum Newton iteration exceeded' or
%             'absolute accuracy reached' or
%             'relative accuracy reached'
%
% This m-file serves as the help file of the maxdet mex-files.

