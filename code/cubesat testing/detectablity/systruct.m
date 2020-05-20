function [nAD]=systruct(W,M,tol)
% SYSTRUCT : System structure determined by two grammians
%
%            [nAD]=systruct(W,M,tol)
%
%            Inputs:
%            W   : Reachability/Controllability grammian
%            M   : Observability/Reconstructability grammian
%            tol : Tolerance (default 1e-12)
%
%            Output:
%            nAD : [nA;nB;nC;nD]
%
% GvW/WdK 30-3-2011         

if nargin==2; tol=1e-12;
elseif nargin~=3; error(' 2 or 3 inputs required'); end
tolmax=1e-2; gk=1; nx=length(W); 

while gk && tol<=tolmax
  % Determine nc
  [U,S,V]=svd(W); tolW = max(size(W)') * max(diag(S)) * tol;
  nc=sum(diag(S)>tolW);

  % Determine no
  [U,S,V]=svd(M); tolM = max(size(M)') * max(diag(S)) * tol;
  no=sum(diag(S)>tolM);

  % Determine nb  
  WM=W*M; [U,S,V]=svd(WM); tolWM = max(size(WM)') * max(diag(S)) * tol;
  nb=sum(diag(S)>tolWM); nb=min([nb,nc,no]);

  % Determine nAD
  nAD=([0 1 0 -1; 0 0 0 1; 1 -1 -1 1; 0 0 1 -1]*[nx;nc;no;nb])';
  if sum(nAD)~=nx || min(nAD)<0 || max(nAD)>nx; tol=10*tol; else gk=0; end
end
if gk; disp(' System structure detection failed'); end