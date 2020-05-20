function [lqgarr,costed,tt,qt,rt,mct,lat,vt,wt]=lqgedat(func,sysdyn,sysout,phi,lqgmatr,dims,satimes,tuxl,xudev,ds,maf,varargin)
% LQGEDAT: Function that computes the equivalent discrete-time LQG problem data
%	         for a synchronous aperiodically sampled linear
%	         time-varying system. All matrices
%	         and references are specified in the function func.
%
%	         [lqgarr,costed,tt,qt,rt,mct,lat,vt,wt]=lqgedat(func,sysdyn,
%           sysout,phi,dims,satimes,tuxl,xudev,ds,maf,varargin)
%
%	Input: 	
%
%	func    : String containing the function name of the function
%	          that  describes the dynamics and the costcriterion
%	          see e.g. aslqgdy.m.
% sysdyn  : String containing function name of state-space model
%           Format [dxdt]=f(t,x,u,varargin) varargin={pv,extinp,...}
% sysout  : String containing function name of state-space output
%           Format [y]=f(t,x,u,varargin) varargin={pv,extinp,...}
% phi     : String containing function name of terminal cost
%           Format [phi]=phi(t,xf,varargin) varargin={pv,extinp,...}
% lqgmatr : Name of file containing the LQG problem data
% dims    :  dimensions
%	satimes : [t0  t1  t2 .. tN;
%	           n0  n1  n2 .. nN]
%	           t0,t1,..,tN: Sampling instances and final time
%	           n0,n1,..,nN: Number of integration steps during sampling period
% tuxl    : [t0  t1 ...  tN;
%            u0* u1* ... uN*;
%            x0* x1* ... xN*
%            l0* l1* ... lN*]' (Transpose!)
%            u0* ... :  optimal control at level 1
%            x0* ... :  optimal state trajectory at level 1
%            l0* ... :  optimal co-state trajectory at level 1
% xudev   :  State and control perturbations to compute numerical derivatives
% ds      :  Scaling factors state
% maf     :  max(abs(f(x,u)))
% varargin:  additional inputs. varargin{1:2}={[pv] [extinp]}
%
%	Output:
%
%	lqgarr:  Equivalent discrete LQG problem data in array format
% costed:  Equivalenr discrete costs
% tt,qt :  t, Q(t)
% rt,mct:  R(t), M(t)
% vt,wt :  V(t), W(t)
% lat   :  lambda(t) from level 1
%
% Array format:
%
%          i=1,2,..,N <-> i'=0,1,..,N-1
%	         ind(i), i=1,2,..,N pointers to starting positions in lqgarr :
%	         the data for time i starts at lqgarr(ind(i)+1) or at
%	         lqgarr(lqgarr(i+1)+1)
%
%          lqgarr=
%	         [N;ind(1);..;ind(N);
%	         n(1);n(2);m(1);l(1);x0m(:);x0c(:);
%	         phi(1)(:);gam(1)(:);q(1)(:);r(1)(:)mc(1)(:);v(1)(:);c(1)(:);w(1)(:);me(1)(:)
%	         n(2);n(3);m(2);l(2);
%	         phi(2)(:);gam(2)(:);q(2)(:);r(2)(:)mc(2)(:);v(2)(:);c(2)(:);w(2)(:);me(2)(:)
%	         ...
%	         ...
%	         n(N-1);n(N);m(N-1);l(N-1);
%	         phi(N-1)(:);gam(N-1)(:);...;v(N-1)(:);
%	         n(N);H(:)];

%
% Initialization
%
dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;

if nargout>=5;
  qt=[]; rt=[]; mct=[]; tt=[]; lat=[]; vt=[]; wt=[];
end;
if nargin==5;tuxl=0;end;
[nl,ml]=size(satimes);N=ml-1;
ni1=0; lqgarr=[]; ind=N+2; indt=ind; costed=0;
hw=waitbar(0,'Equivalent discrete-time LQG problem computation');
[x0m,x0c,H]=lqgx0gh(phi,lqgmatr,dims,tuxl(end,1),tuxl(end,2+dimu:1+dimu+dimxd)',xudev,ds);

%
%	Computation of recursive solution
%
flag=0; save 'lqgedat.txt' flag -ascii;
for i=1:N
  tis=satimes(1,i);tif=satimes(1,i+1);ni=satimes(2,i);
  ts=tif-tis;time=[tis:(tif-tis)/ni:tif];
  u=tuxl(i,2:dimu+1)'; x=tuxl(i,2+dimu:1+dimu+dimxd)';
%
%	Obtain continuous dynamics
%
  [ah,bh,vh,qh,mch,rh,xrh,urh,dh,n,m,c,w,l,la]=feval(func,sysdyn,sysout,lqgmatr,dims,time,tuxl,xudev,ds,maf,varargin{:});
	if nargout>=3;
    qt=[qt qh]; rt=[rt rh]; mct=[mct mch]; tt=[tt time]; lat=[lat la];
    vt=[vt vh]; wt=[wt w(:)];
  end;
%
%	Convert to equivalent discrete-time problem and store
%
  [phi,gam,qi,ri,mc,vi,gk]=edortv(ah,bh,qh,rh,vh,ts,mch); n=size(phi,1);
  me=zeros(n,l); costed=costed+gk;

% Symmetrize and ensure nonnegativeness
  [qi,mc,ri]=pm2snn(qi,mc,ri); [vi]=cin2nn(vi);

% Adjust qi based on second order terms output equation
  inx=[1:dimxd+dimu+1];
  [q1]=slqmr(sysout,dims,tis,[x;0],u,xudev(inx),ds,'g',1,varargin{:});
  qi(1:dimxd,1:dimxd)=qi(1:dimxd,1:dimxd)+ts*q1(1:dimxd,1:dimxd);
  
  [ni,mi,nh,ldi]=pgcvchk(phi,gam,qi,ri,mc,vi,c,w,me);
  if i~=1;
    if ni~=ni1; error('  ni sequence error'); end;
    ind=ind+4+ni*nh+nh*mi+ni*ni+mi*mi+ni*mi+nh*nh+ldi*ni+ldi*ldi+nh*ldi;
    ni1=nh; lqgarr=[lqgarr; ni; ni1; mi; ldi; phi(:); gam(:); qi(:); ri(:); mc(:); vi(:); c(:); w(:); me(:)];
  else
    [nn,mm]=size(x0m);
    if mm~=1; error(' x0m should be a column vector'); end
    if nn~=ni error(' Dimensions x0m incompatible with p(1)'); end
    [nn,mm]=size(x0c);
    if nn~=dimxd || mm~=dimxd; error(' Dimensions x0c incompatible a'); end
    x0c=[x0c zeros(n,ni-n); zeros(ni-n,ni)];
    ind=ind+4+ni+ni*ni+ni*nh+nh*mi+ni*ni+mi*mi+ni*mi+nh*nh+ldi*ni+ldi*ldi+nh*ldi;
    ni1=nh; lqgarr=[lqgarr; ni; ni1; mi; ldi; x0m(:); x0c(:); phi(:); gam(:); qi(:); ri(:); mc(:); vi(:); c(:); w(:); me(:)];
  end
  indt=[indt; ind]; waitbar(i/N);
  if fileflag('lqgedat.txt')
     close(hw); error(' lqgedat.txt invoked termination');
  end;
end;
close(hw);
[nn,mm]=size(H);
if nn~=ni1 || mm~=ni1;
   error(' Dimensions H incompatible with p(N-1)');
end;
lqgarr=[N+1;indt;lqgarr;ni1;H(:)];