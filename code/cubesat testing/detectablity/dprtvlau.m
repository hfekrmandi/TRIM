% DPRTVLAU:Deterministic Parameter Reduced Order Compensation.
%          Optimal control of (P(i),G(i),C(i),V(i),W(i),x0m,x0c)
%          based on UDU factored Lyapunov equations
%          for asynchronously sampled systems
%          (Q(i),R(i),H) and nc<=nx.
%          (P(i),C(i)) deterministic.
%          Time shift i=0,1,..,N -> i=1,2,..,N+1
%
%          [xfkl,sigp,sigs,ppt,sst,nsweeps]=dprtvlau(lqgarr,nc,opt);
%          or
%          [xfkl,sigp,sigs,ppt,sst,nsweeps]=dprtvlau(lqgarr,nc,opt,ppt,sst);
%
% Input:
%          lqgarr: column vector containing the lqg problem data
%          nc    : row vector with dimensions of the compensator state
%          opt   : optional array containing algorithm parameter settings.
%          opt(1): convergence tolerance, epsl, default 1e-6. 
%          opt(2): iteration limit, maxtps, default 1e20.
%          opt(3): maximum number of iterations (sweeps), itm, default 5000. 
%          opt(4): plot parameter, itc, default 0.
%                  0    : no plot.
%                  <0   : plot at the end
%                  >0   : plot every itc iterations.
%          opt(5): numerical damping coefficient, 0<b<1, default 0.25.
%                  if opt(5)>1 | opt(5)<0 then b=0
%          opt(6): initial conditions ph, sh, ic=1,2,3, default 1
%                  ic=1: 0.9*eye(nx)+0.1*ones(nx)
%                  ic=2: random
%                  ic=3: eye(nx)
%          ppt,sst: optional if iteration is to be continued.
%
% Output:
%          Time shift i=0,1,..,N -> i=1,2,..,N+1
%          xfkl : Optimal compensator in array format
%          sigp,sigs: Minimum costs (if sigp=sigs)
%          ppt, sst: Solution lyapunov equations.
%          nsweeps : number of forward and backward iterations
%
% Comments:
%          See also dproeqlf dproeqlb gmhfac
%
%          L.G. Van Willigenburg, W.L. De Koning, 27-06-96.
%
  function [xfkl,sigp,sigs,ppt,sst,nsweeps]=dprtvlau(lqgarr,nc,opt,ppt,sst);

%
% Check number of input arguments and get
% problem data from lqgarr
%
  if nargin==2; no=0;
  elseif nargin==3; no=max(size(opt));
  elseif nargin==5; no=max(size(opt));
  else; error('two,three or five input arguments required'); end;
%
% Check data
%
  [nh,mh]=size(lqgarr); N=lqgarr(1);
  if mh~=1 ;
    error(' First input argument must be a column vector');
  end;
  [nh,mh]=size(nc);
  if nh~=1 | mh<N-1 ;
    error(' Second input argument must be a row vector with length >= N-1');
  end; nc(N)=0; nc=nc(1:N);
%
% Algorithm parameter settings.
%
  opth=zeros(1,6);
  if no~=0; opth(1:no)=opt(1:no); end;
  opt=opth;
  if opt(1)>0; epsl=opt(1); else; epsl=1e-6; end
  if opt(2)>0; maxtps=opt(2); else; maxtps=1e20; end
  if opt(3)>0; itm=opt(3); else; itm=5000; end
  itc=opt(4);
  if opt(5)==0; b=0.25;
  elseif opt(5)<0 | opt(5)>1;
  b=0; else b=opt(5); end; bb=b;
  if opt(6)>=1 & opt(6)<=3; ic=round(opt(6)); else; ic=1; end
  
% Dimension check and index array computation
  nxt=0; inps=nxt;
  for i=1:N
    [nx]=getedo(lqgarr,i); nc(i)=min(nc(i),nx);
    nxt=nxt+(nx+nc(i))*(nx+nc(i)); inps=[inps; nxt];
    if nc(i)<0 | nc(i) > nx;
       error([' 0<=nc(' num2str(i) ')<=n(' num2str(i) ') violated']);
    end;
  end
  
% Boundary conditions
  [x0m,x0c,H]=getxgh(lqgarr);

% Initialization PP, SS if necessary
  if nargin~=5
    [ppt,sst]=inipstva(lqgarr,nc,ic,1);
  end;

% Initialization other variables
  nsweeps=0; endloop=0; conver=0; tr=0; trr=0; trt=[];

%
% SS, PP loop
%
  while endloop~=1

    nsweeps=nsweeps+1; tro=tr;
    if bb~=0; ppto=ppt; ssto=sst; end;

%   SS loop
    sigst=0;
    pn1=inps(N)+1; pn2=inps(N+1);
    [nx]=getedo(lqgarr,N);
    ss=zeros(nx+nc(N)); ss(:)=sst(pn1:pn2,:);
    for i=N-1:-1:1
      [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me]=getedo(lqgarr,i);
      pn1=inps(i)+1; pn2=inps(i+1);
      [nx]=getedo(lqgarr,i);
      pp=zeros(nx+nc(i)); pp(:)=ppt(pn1:pn2,:);
      [pp,ss,f,k,l,sigp,sigs,tps]=dprludub(pp,ss,p,g,c,q,r,v,w,mc,me,bb);
      sst(pn1:pn2,:)=ss(:); sigst=sigst+sigs;
    end
    sigst=sigst+trace(ud2sym(ss)*ud2sym(pp));

%   Numerical damping
    if bb~=0;
      for i=N-1:-1:1
        [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me]=getedo(lqgarr,i);
        pn1=inps(i)+1; pn2=inps(i+1);
        [nx]=getedo(lqgarr,i);
        ss=zeros(nx+nc(i)) ; ss(:)=sst(pn1:pn2,:);
        sso=zeros(nx+nc(i)); sso(:)=ssto(pn1:pn2,:);
        [ss]=uduadd(ss,sso,1-bb,bb,1e-12,-1e12);
        sst(pn1:pn2,:)=ss(:);
      end
    end

    [nx]=getedo(lqgarr,1); nxc=nx+nc(1);
    s12=ss(1:nx,nx+1:nxc); s2=ss(nx+1:nxc,nx+1:nxc);
    if nxc~=nx
      [uh]=udinv(s2,1e-12,-1e12); [uh]=ud2ud(uh); hh=-uh'*s12';
    else
      hh=zeros(0,nx);
    end
    x0e=hh*x0m;
    p1=x0m*x0m'+x0c; p12=x0m*x0e'; p2=x0e*x0e';
    pp=[p1 p12; p12' p2]; pp=sym2ud(pp,1e-12,-1e12);
    ppt(pn1:pn2,:)=pp(:);

%   PP loop
    sigpt=0;
    for i=1:N-1
      [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me]=getedo(lqgarr,i);
      pn1=inps(i+1)+1; pn2=inps(i+1+1);
      ss=zeros(nxn+nc(i+1)); ss(:)=sst(pn1:pn2,:);
      [pp,ss,f,k,l,sigp,sigs,tps]=dprluduf(pp,ss,p,g,c,q,r,v,w,mc,me,bb);
      ppt(pn1:pn2,:)=pp(:); sigpt=sigpt+sigp;
    end
    sigpt=sigpt+trace(ud2sym(ss)*ud2sym(pp));

%   Numerical damping
    if bb~=0;
      for i=1:N
        [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me]=getedo(lqgarr,i);
        pn1=inps(i)+1; pn2=inps(i+1);
        [nx]=getedo(lqgarr,i);
        pp=zeros(nx+nc(i)) ; pp(:)=ppt(pn1:pn2,:);
        ppo=zeros(nx+nc(i)); ppo(:)=ppto(pn1:pn2,:);
        [pp]=uduadd(pp,ppo,1-bb,bb,1e-12,-1e12);
        ppt(pn1:pn2,:)=pp(:);
      end
    end

%   Convergence computations
    pn1=inps(1)+1; pn2=inps(2);
    [nx]=getedo(lqgarr,1); nxc=nx+nc(1);
    ss=zeros(nxc); ss(:)=sst(pn1:pn2,:);
    [st,sh]=ud2psh(ss,nc(1));
    tr=trace(st);
%    tr=trace(ud2sym(ss(1:nx,1:nx)));

    pn1=inps(N)+1; pn2=inps(N+1);
    [nx]=getedo(lqgarr,N); nxc=nx+nc(N);
    pp=zeros(nxc); pp(:)=ppt(pn1:pn2,:);
    [pt,ph]=ud2psh(pp,nc(N));
    tr=tr+trace(pt);
%    tr=tr+trace(ud2sym(pp(1:nx,1:nx)));

    trt=[trt tr]; trr=abs((tr-tro)/(tr+tro));
%
% Test convergence.
%
    if nsweeps>itm
      disp('Maximum number of iterations (sweeps) exceeded');
      endloop=1; sigp=Inf; sigs=sigp; sigc=sigp;
    elseif trr>maxtps
      disp('Algorithm failed, divergence');
      endloop=1; sigp=Inf; sigs=sigp; sigc=sigp;
    elseif trr<epsl;
      conver=conver+1;
    end

    if conver==2; bb=0.9*b;
    elseif conver==3; endloop=1; sigp=0;
    else; bb=b; end;

%    if nsweeps==100; endloop=2; end;
%
% Plot progress computation.
%
    if (itc>0 & rem(nsweeps,itc)==0) | (itc~=0 & endloop==1)
      plconv(trt,trr,epsl,endloop,'DPRTVLAU','Trace(P1(N)+S1(0))',nx,nc);
    end;

    [flag]=fileflag('stop.txt');
    if flag; delete stop.txt; error(' Termination due to stop.txt'); end
    
  end;
  delete stop.txt

  if ~isinf(sigp)
%
% Final update 
%
    ind=N+1; indt=ind; xfkl=[];
    pn1=inps(1)+1; pn2=inps(1+1);
    [nx]=getedo(lqgarr,1); nxc=nx+nc(1);

    ss=zeros(nxc); ss(:)=sst(pn1:pn2,:);
    s12=ss(1:nx,nx+1:nxc); s2=ss(nx+1:nxc,nx+1:nxc);
    if nxc~=nx
      [uh]=udinv(s2,1e-12,-1e12); [uh]=ud2ud(uh); hh=-uh'*s12';
    else
      hh=zeros(0,nx);
    end
    x0e=hh*x0m;

    pp=ss; pp(:)=ppt(pn1:pn2,:);
    for i=1:N-1
      [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me]=getedo(lqgarr,i);
      pn1=inps(i+1)+1; pn2=inps(i+1+1);
      ss=zeros(nxn+nc(i+1)); ss(:)=sst(pn1:pn2,:);
      [pp,ss,f,k,l]=dprluduf(pp,ss,p,g,c,q,r,v,w,mc,me,0);
      if i~=1
        xfkl=[xfkl; nc(i); nc(i+1); mi; ldi; f(:); k(:); l(:)];
        ind=ind+4+nc(i)*nc(i+1)+nc(i+1)*ldi+nc(i)*mi;
      else
        xfkl=[xfkl; nc(i); nc(i+1); mi; ldi; x0e(:); f(:); k(:); l(:)];
        ind=ind+4+nc(i)+nc(i)*nc(i+1)+nc(i+1)*ldi+nc(i)*mi;
      end
      indt=[indt; ind];
    end
    xfkl=[N;indt;xfkl;nc(N)];
    sigp=sigpt; sigs=sigst;
  else;
    sigc=Inf; sigs=sigp; xfkl=[];
  end
