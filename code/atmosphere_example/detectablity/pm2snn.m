function [qi,mc,ri]=pm2snn(qi,mc,ri);
% PM2SNN : Function to symmetrize and make nonnegative the partioned matrix [qi mc; mc' ri]
%
%         Input :
%                qi,mc,ri
%
%         Output:
%                qi,mc,ri
%
% GvW 10-5-2005

if nargin~=3; error(' 3 inputs required'); end;
if nargout~=3; error(' 3 outputs required'); end;

[n,n1]=size(qi); [m,m1]=size(ri); [n2,m2]=size(mc);

if n1~=n; error(' qi must be square'); end;
if m1~=m; error(' ri must be square'); end;
if n2~=n | m2~=m; error(' mc incompatible with qi, ri'); end;

qs=cin2nn([qi mc; mc' ri]);

qi=qs(1:n,1:n); mc=qs(1:n,n+1:n+m); ri=qs(n+1:n+m,n+1:n+m);
