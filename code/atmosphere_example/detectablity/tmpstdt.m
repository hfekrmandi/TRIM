function [npt,nst,ptt,stt,NW,NM,rW,rM,nAD,hfs]=tmpstdt(lqgarr,js,ks,tolWM,epsr,epsw)
% TMPSTDT : Temporal and one-step stabilizability and detectability analysis
%           of a time-varying linear(ized) discrete-time system
%           Described in paper:
%           L.G. Van Willigenburg, W.L. De Koning, 2013,
%           "Temporal and one-step stabilisability and detectability of
%            discrete time linear systems",
%           IET Control Theory & Applications, 7, 1, 151-159
%           Download this paper from:
%           http://gvw007.yolasite.com/resources/TempOneIET2013.pdf
%
%           [npt,nst,ptt,stt,NW,NM,rW,rM,nAD,hfs]=tmpstdt(lqgarr,js,ks,tolWM,epsr,epsw)
%
%           Inputs:
%           lqgarr  : Array representing a discrete-time LQG problem
%                     from which the system matrices are extracted
%                     for analysis
%           js      : J-steps for J-step controll./reachability analysis
%           ks      : K-steps for K-step reconstruct./observability analysis
%           tolWM   : relative tolerance to determine rank(W), rank(M)
%                     W: J-step controll./reachability grammian
%                     M: K-step reconstruct./observability grammian
%           epsr    : R=epsrw*eye(nu); R value temporal stabilizability analysis
%           epsw    : W=epsrw*eye(ny); W value temporal detectability analysis
%
%           Outputs:
%           npt     : Norms pt over temporal unreachable/uncontrollable int
%           nst     : Norms st over temporal unreconstr./unobservable int
%           ptt     : Values pt over temporal unreachable/uncontrollable int
%           stt     : Values st over temporal unreconstr./unobservable int
%           NW      : J-step unreachable times
%           NM      : K-step unobservable times
%           rW      : Rank(W(i)), possibly computed from nAD, see below
%           rM      : Rank(M(i)), possibly computed from nAD, see below
%           nAD     : Temporal linear system structure from
%                     differential Kalman decomposition (may fail).
%                     If not specified rank(W(i)) and rank(M(i))
%                     are computed through rankr which does not fail.
%           hfs     : handles figures
%
% Gerard van Willigenburg / Willem de Koning 20-4-2011

% Check inputs
[ni,ni1,mi,ldi,phi,gam,qi,ri,mc,vi,c,w,me,x0m,x0c]=getedo(lqgarr,1);
N=lqgarr(1);
tolWMd=1e-10; epsrwd=1e-10; %Default tolerances

if nargin==1; js=2*ni; ks=2*ni; tolWM=tolWMd; epsr=epsrwd; epsw=epsrwd;
elseif nargin==3; tolWM=tolWMd; epsr=epsrwd; epsw=epsrwd;
elseif nargin==4; epsr=epsrwd; epsw=epsrwd;
elseif nargin~=6 error(' Wrong number of inputs'); end 

% Get system matrices from lqgarr
for k=1:N
  [ni,ni1,mi,ldi,phi,gam,qi,ri,mc,vi,c,w,me,x0m,x0c]=getedo(lqgarr,k);
  Fo{k}=phi; Go{k}=gam; Ho{k}=c;
end

%% Determine temporal linear j-step k-step structure

% Determine grammians
nx=length(Fo{1}); Wn{1}=zeros(nx); Mn{N}=zeros(nx);
for k=1:N-1
% Grammians over [k-nx+1,k+nx-1]
  Wn{k+1}=zeros(nx); Mn{k}=zeros(nx);
  for j=max(1,k-js+1):k
    Wn{k+1}=Fo{j}*Wn{k+1}*Fo{j}'+Go{j}*Go{j}'; Wn{k+1}=0.5*(Wn{k+1}+Wn{k+1}'); 
    %Fd=Fo{j}, Gd=Go{j}, Wd=Wn{k+1}, svd(Wd)
  end
  for j=min(N-1,k+ks-1):-1:k;
    Mn{k}=Fo{j}'*Mn{k}*Fo{j}+Ho{j}'*Ho{j}; Mn{k}=0.5*(Mn{k}+Mn{k}');
  end
end

%% Determine temporal linear system structure
nAD=zeros(N,4); rW=zeros(N,1); rM=rW;
for k=1:N;
  if nargout>8;
    nAD(k,:)=systruct(Wn{k},Mn{k},tolWM);
    rW(k)=sum(nAD(k,1:2)); rM(k)=sum(nAD(k,[2 4]));
  else
    rW(k)=rankr(Wn{k},tolWM); rM(k)=rankr(Mn{k},tolWM);
  end
end
if nargout>8; disp('     i    nAD'); disp([[0:N-1]' nAD]);
else disp('     i Rank(W) Rank(M)'); disp([[0:N-1]' rW rM]); end

%% Compute temporal and differential stabilizability
%% and detectability measures from the two associated LQ problems

%%Compute time instants where the system is not js reachable
[m,NW]=stackomp(rW(1+js:end),'<',nx,1); NW=NW+js;

% Test values split
%NW=[1:4 6:14 16 19:25 29 33:36 38:40]';

% Split NW into intervals
NWn={}; nW=0; NWr=length(NW); i1=1; i2=i1;
while i2<NWr
  while i2<NWr && NW(i2)==NW(i2+1)-1; i2=i2+1; end
  if i2>=i1+js; nW=nW+1; NWn{nW}=NW(i1:i2-js); end
  i1=i2+1; i2=i1;
end

% Initialization
q=zeros(ni); r=epsr*eye(mi); mc=zeros(ni,mi);
stt=[]; nst=[]; nsi=[]; NW=[]; NWi=[];

for k=nW:-1:1
  % LQ problem associated with temporal and differential stabilizability
  NWt=NWn{k};
  if length(NWt)>1
    NW=[NWt; NW]; NWi=[NWt(end); NWi];
    st=eye(ni); stn=norm(st); nst=[stn nst]; nsi=[stn nsi];
    for k=fliplr(NWt(2:end)')
      [kt,st]=sric1(Fo{k},Go{k},q,r,mc,st);
      stt=[st stt]; stn=norm(st); nst=[stn nst];
    end
  end
end

%%Compute discrete time instants where the system is not ks observable
[m,NM]=stackomp(rM(1:end-ks),'<',nx,1);
 
% Test values split
%NM=[1:4 6:14 16 19:25 29 33:36 38:40]';

% Split NM into intervals
NMn={}; nM=0; NMr=length(NM); i1=1; i2=i1;
while i2<NMr
  while i2<NMr && NM(i2)==NM(i2+1)-1; i2=i2+1; end
  if i2>=i1+ks; nM=nM+1; NMn{nM}=NM(i1+ks:i2); end
  i1=i2+1; i2=i1;
end

% Initialization
v=zeros(ni); w=epsw*eye(mi); me=zeros(ni,ldi);
ptt=[]; npt=[]; npi=[]; NM=[]; NMi=[];

for k=1:nM
  % LQ problem associated with temporal and differential stabilizability
  NMt=NMn{k};
  if length(NMt)>1
    NM=[NM; NMt]; NMi=[NMi; NMt(1)];
    pt=eye(ni); ptn=norm(pt); npt=[npt ptn]; npi=[npi ptn];
    for k=NMt(2:end)'
     [lt,pt]=sric1(Fo{k}',Ho{k}',v,w,me,pt);
     ptt=[ptt pt]; ptn=norm(pt); npt=[npt ptn];
    end
  end
end

%% Plot
jstr=num2str(js); kstr=num2str(ks);
if ~isempty(NW);
  hfs=figure; plot(NW'-1,nst,'o'); legend([jstr '-step uncontroll./unreach. times'])
  title([' Temporal and one-step stabilizability measures']);
  xlabel('Discrete-time'); ylabel('||S_{i}^{\epsilon}||');
  % Plot initial and final points intervals
  if ~isempty(NWi); hold on; plot(NWi-1,nsi,'*r'); hold off; end;
else
  disp([' No instants where the linear(ized) system is ' jstr '-step uncontroll./unreach.']);
end
if ~isempty(NM);
  hf=figure; hfs=[hfs hf]; plot(NM'-1,npt,'o'); legend([kstr '-step unreconstr./unobserv. times'])
  title([' Temporal and one-step detectability measures']);
  xlabel('Discrete-time'); ylabel('||P{i}^{\epsilon}||');
  % Plot initial and final points intervals
  if ~isempty(NMi); hold on; plot(NMi-1,npi,'*r'); hold off; end;
else
  disp([' No instants where the linear(ized) system is ' kstr '-step unreconstr./unobserv.']);
end
