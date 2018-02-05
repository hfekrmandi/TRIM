% GETEDO : Get equivalent discrete-time system matrices
%          for time i from the edorp array lqgarr.
%
%          [ni,ni1,mi,ldi,phi,gam,qi,ri,mc,vi,c,w,me,x0m,x0c]=getedo(lqgarr,i)
%
%	   input:
%          lqgarr : edorp array (see lqgrarv, lqgdat)
%          i      : time index 1,2,..,N (N=N+1)
%
%          output :
%          ni,ni1,mi,ldi
%          phi,gam,qi,ri,am,vi at time i
% 
%    Comments:
%          See also: lqgdat, lqgrarv
%
%          L.G. Van Willigenburg, W.L. De Koning, 28-11-95.
%
function [ni,ni1,mi,ldi,phi,gam,qi,ri,mc,vi,c,w,me,x0m,x0c]=getedo(lqgarr,i);

if nargin~=2; error('  Two input arguments required'); end;
[nh,mh]=size(lqgarr);
if mh~=1; error(' First input argument must be a column vector'); end;

% Determine N (=N+1) and check i
N=lqgarr(1); if i<1 | i>N; error(' 1 <= i <= N violated'); end;

% Adjust i to N-1 if i==N :
% then only n(N) is important
if i==N; hN=1; i=N-1; else; hN=0; end;

% Determine the offset o(i)
oi=lqgarr(1+i);

% Determine n(i),n(i+1),m(i) and
oi=oi+1; ni=lqgarr(oi); oi=oi+1; ni1=lqgarr(oi);
oi=oi+1; mi=lqgarr(oi); oi=oi+1; ldi=lqgarr(oi);

% Read x0m,x0c if i==1
if i==1;
	x0m=zeros(ni,1); x0m(:)=lqgarr(oi+1:oi+ni); oi=oi+ni;
	x0c=zeros(ni); x0c(:)=lqgarr(oi+1:oi+ni*ni); oi=oi+ni*ni;
else
	x0m=[]; x0c=[];
end

% Obtain matrices
n=ni1; m=ni; nm=n*m;
phi=zeros(n,m); phi(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=mi; nm=n*m;
gam=zeros(n,m); gam(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ni; m=ni; nm=n*m;
qi=zeros(n,m); qi(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=mi; m=mi; nm=n*m;
ri=zeros(n,m); ri(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ni; m=mi; nm=n*m;
mc=zeros(n,m); mc(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=ni1; nm=n*m;
vi=zeros(n,m); vi(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ldi; m=ni; nm=n*m;
c=zeros(n,m); c(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ldi; m=ldi; nm=n*m;
w=zeros(n,m); w(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=ldi; nm=n*m;
me=zeros(n,m); me(:)=lqgarr(oi+1:oi+nm); oi=oi+nm;
if lqgarr(oi+1)~=ni1; error(' ni1 check error'); end;

% Adjust ni(N) if i was equal to N 
if hN==1; ni=ni1; end;
