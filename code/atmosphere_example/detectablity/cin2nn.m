function [snn]=cin2nn(s,mt)
% CIN2NN: Convert square indefinite matrix to a square non-negative matrix
%         Application: computation of Q,R,M in LQG problems from d2H/dx, d2H/du
%         or d2f/dx2, d2f/du2
%
%         [snn]=cin2nn(s,m,tol)
%
%         s   : square indefinite input matrix
%         snn : square non-negative symmetric output matrix
%         mt  : method 1: Nonnegative definite upperbound from SVD
%                      2: Nonnegative definite diagonal upperbound
%                      3: Nonnegative definite diagonal matrix, but not an upperbound
%         snn : square non negative definite matrix
%         

if nargin<2; mt=1; tol=1e-6;
elseif nargin<3; tol=1e-6;
elseif tol<=0; error(' tol <= 0'); end;

% Check squareness and symmetry of s
[n,m]=size(s);

if n~=m;
  error(' s must be square');
end;

% Symmetrize s
s=0.5*(s+s');

if mt==1
% Nonnegative symmetric upperbound from SVD
  [u,snn,v]=svd(s); snn=u*snn*u';
   
elseif mt==2
% Nonnegative diagonal upperbound
  snn=cumsum(abs(s)); snn=diag(snn(end,:));
      
elseif mt==3
% Diagonal matrix, no upperbound
  snn=max(abs(s)); snn=diag(snn(end,:));
else
  error('  Second input must be 1,2 or 3');
end