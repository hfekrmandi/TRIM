% GETXGH : Get initial state and final penalty data
%
%          [x0m,x0c,H]=getxgh(pgqrav);
%
%	   input:
%          pgqrav : edorp array (see lqgrarv)
%
%          output : x0m,x0c,H
% 
function [x0m,x0c,H]=getxgh(pgqrav);

if nargin~=1; error('  One input argument required'); end;
[nh,mh]=size(pgqrav);
if mh~=1; error(' First input argument must be a column vector'); end;

% Determine N and the offset o(i) to obtain x0m,x0c
N=pgqrav(1); oi=pgqrav(1+1);

% Determine n(i),n(i+1),m(i) and
oi=oi+1; ni=pgqrav(oi); oi=oi+1; ni1=pgqrav(oi);
oi=oi+1; mi=pgqrav(oi); oi=oi+1; ldi=pgqrav(oi);

% Obtain matrices
n=ni; m=1; nm=n*m;
x0m=zeros(n,m); x0m(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni; m=ni; nm=n*m;
x0c=zeros(n,m); x0c(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=ni; nm=n*m;

% Get remaining matrices for i==1 to check ni1
phi=zeros(n,m); phi(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=mi; nm=n*m;
gam=zeros(n,m); gam(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni; m=ni; nm=n*m;
qi=zeros(n,m); qi(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=mi; m=mi; nm=n*m;
ri=zeros(n,m); ri(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni; m=mi; nm=n*m;
mc=zeros(n,m); mc(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=ni1; nm=n*m;
vi=zeros(n,m); vi(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ldi; m=ni; nm=n*m;
c=zeros(n,m); c(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ldi; m=ldi; nm=n*m;
w=zeros(n,m); w(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
n=ni1; m=ldi; nm=n*m;
me=zeros(n,m); me(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
if pgqrav(oi+1)~=ni1; error(' ni1 check error'); end;

% Determine the offset o(i) to obtain H
oi=pgqrav(1+N);

% Determine n(i),n(i+1),m(i) and
oi=oi+1; ni=pgqrav(oi);
n=ni; m=ni; nm=n*m;
H=zeros(n,m); H(:)=pgqrav(oi+1:oi+nm); oi=oi+nm;
if oi~=max(size(pgqrav)); error(' H check error'); end;
