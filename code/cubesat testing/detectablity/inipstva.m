% INIPSTVA: Initialization finite-horizon
%           Lyapunov equations and SDOPE
%           with time-varying dimensions
%
%           [ppt,sst,ptt,stt,pht,sht]=inipstva(lqgarr,nc,ic,udu);
%           [ptt,stt,pht,sht]=inipstva(lqgarr,ic,udu);
%           [ppt,sst]=inipstva(lqgarr,ic,udu);
%
%           Input:
%           lqgarr : LQG problem data in array format
%           nc     : Prescribed dimensions compenmsator
%           ic     : flag for type of initialization
%           udu    : optional if 1 then UDU factored PP, SS
%
%           Output:
%           ppt, sst : PP(1..N), SS(1..N)
%           ptt, stt : P(1..N), S(1..N)
%           pht, sht : Phat(1..N), Shat(1..N)
%
%           GvW/WdK 22-2-02
%
function [ppt,sst,ptt,stt,pht,sht]=inipstva(lqgarr,nc,ic,udu);

if nargin<3; error(' 3 or 4 input arguments required');
elseif nargin==3; udu=0; end;
N=lqgarr(1); i=1; [nx]=getedo(lqgarr,i);

[x0m,x0c,H]=getxgh(lqgarr);
x0e=ones(nc(1),1);
p1=x0m*x0m'+x0c; p12=x0m*x0e'; p2=x0e*x0e';
ppb=[p1 p12; p12' p2]; [ptb,phb]=ps2psh(ppb,nx,nc(i));
if udu; ppb=sym2ud(ppb,1e-12,-1e12); end;

[nx]=getedo(lqgarr,N);
zxc=zeros(nx,nc(N)); zc=zeros(nc(N));
ssf=[H zxc; zxc' zc]; [stf,shf]=ps2psh(ssf,nx,nc(N));
if udu; ssf=sym2ud(ssf,1e-12,-1e12); end;

ppt=ppb(:); ptt=ptb(:); pht=phb(:);
for i=2:N
  [nx]=getedo(lqgarr,i);
  [ph,p12,p2]=iniphsh(ic,nx,nc(i));
  pp=[ph p12; p12' p2]; [pt,ph]=ps2psh(pp,nx,nc(i));
  if udu; pp=sym2ud(pp,1e-12,-1e12); end
  ppt=[ppt; pp(:)]; ptt=[ptt; pt(:)]; pht=[pht; ph(:)];
end
sst=[]; stt=[]; sht=[];
for i=1:N-1
  [nx]=getedo(lqgarr,i);
  [sh,s12,s2]=iniphsh(ic,nx,nc(i));
  ss=[sh s12; s12' s2]; [st,sh]=ps2psh(ss,nx,nc(i));
  if udu; ss=sym2ud(ss,1e-12,-1e12); end;
  sst=[sst; ss(:)]; stt=[stt; st(:)]; sht=[sht; sh(:)];
end
sst=[sst; ssf(:)]; stt=[stt; stf(:)]; sht=[sht; shf(:)];

% Adaptation for call with 4 arguments
if nargout==4; ppt=ptt; sst=stt; ptt=pht; stt=sht; end;