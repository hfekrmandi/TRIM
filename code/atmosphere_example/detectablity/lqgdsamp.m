function [a,b,v,q,mc,r,xr,ur,dimd,dimx,dimu,c,w,dimy,la]=lqgdsamp(sysdyn,sysout,lqgmatr,dims,tim,tuxl,xudev,ds,maf,varargin)
%
% LQGDSAMP : Function generates the Digital LQG compensator data for each sampling interval
%            based on the linearized dynamics about the optimal trajectory specified by tuxl.
%
%            L.G. van Willigenburg 31-1-'96.
%
%            [a,b,q,r,xr,ur,d,n,m,c,w,l,h,G]=lqgdsamp(sysdyn,sysout,dims,tim,tuxl)
%
%	Function generates matrix arrays a,v,b,c,w,q,r,xr,ur,d
%	which belong to a digital LQG problem for the time vector
%	tim. The matrix arrays are stored columnwise! See getmat.m!!
%
% Input
%     sysdyn : String with name of function containing the system dynamics
%     sysout : String with name of function containing the system output
%    lqgmatr :  Name of file containing the LQG problem data
%     dims   : dimensions
%     tim    : time vector
%     tuxl   : [t0  t1 ...  tN;
%               u0* u1* ... uN*;
%               x0* x1* ... xN*
%               l0* l1* ... lN*]' (Transpose!)
%               u0* ... :  optimal control at level 1
%               x0* ... :  optimal state trajectory at level 1
%               l0* ... :  optimal co-state trajectory at level 1
%     xudev  : state & control perturbations for finte differences
%     ds     : scaling factors
%     maf    : max(|f|)
% varargin   : additional inputs. varargin{1:2}={[pv] [extinp]}
%
%	Output:
%
%	    a: Matrix array stored columnwise.
%	    the same for b,v,q,r,xr,ur,d. See getmat.m.
%
%	    dimx :	dimension state vector
%	    dimu :	dimension control vector
%
%	    c    :	output matrix at initial time; c(tis)
%	    w    :	covariance matrix of measurement noise w(tis)
%	    dimy :	dimension of the output
%
%	    h    :	matrix representing the quadratic cost at the final time
%	    G    :	matrix representing the initial state covariance
%

  dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;

  if nargout<15 || nargout>16 error('lqgdsamp should have 15 to 16 output arguments'); end;
  if nargin<6; error('lqgdsamp should have 6 or more input arguments'); end;
  [nt,mt]=size(tim); if nt~=1;error('argument should be a row vector'); end;
  
  a=[];b=[];q=[];r=[];mc=[];xr=[];ur=[];d=[];v=[]; la=[];

% Determination initial (tis) and final (tif) time sampling period,
% sampling interval (ts) and stepsize integration (dt)
  tis=tim(1);tif=tim(end);
  [t,ni]=size(tim); ni=ni-1; dt=(tif-tis)/ni;
  [t,ind]=min(abs(tuxl(:,1)-tis));
  x0=[tuxl(ind,2+dimu:1+dimu+dimxd)';0];
  la=[tuxl(ind+1,2+dimu+dimxd:1+dimu+dimxd+dimxd)';1];

  if mt~=1

%   Integration from tis to tif to determine X(t)
    tk=[tis; tif]; uk=tuxl(ind,2:1+dimu);
    [xk,yk,t,x]=digsim(sysdyn,dims,tk,uk,x0,dt,sysout,0,varargin{:});

%   Determination linearized model A(t),B(t) about X(t),U(t) and
%   Q(t),R(t),V(t),D(t),XR(t),UR(t) the other LQG problem data
    na=1e6; 
    for i=1:ni+1
      [ah,bh,qh,rh,mh,vh]=lqgabqrv(sysdyn,lqgmatr,dims,t(i),x(i,:)',tuxl(ind,2:1+dimu)',xudev,ds,la,maf,varargin{:});
      a=[a ah(:)]; b=[b bh(:)]; mc=[mc mh(:)];
      q=[q qh(:)]; r=[r rh(:)]; v=[v vh(:)];
    end
    nh=size(ah,1); xr=zeros(nh,ni+1); ur=zeros(dimu,ni+1);
    
  end

% Generation C,W(tis),x0m,H,G
 [c,w]=lqgcw(sysout,lqgmatr,dims,tis,x(1,:)',tuxl(ind,2:1+dimu)',xudev,ds,la,varargin{:});