function [A,B,Q,R,M,V]=lqgabqrv(sysdyn,lqgmatr,dims,t,x,u,xudev,ds,la,maf,varargin)
% LQGABQRV: LQG compensator design matrices A,B,Q,R,V
%            
%           [A,B,Q,R,V]=lqgabqrv(t,x,u,xpert,upert)
%
%           Input:
%                sysdyn :  string containing filename state=space model
%                 dims  :  dimensions
%               lqgmatr :  Name of file containing the LQG problem data
%                 t     :  current time
%                 x     :  current state including running costs (Mayer formulation)
%                 xudev :  control and state perturbations for derivative computation
%                 ds    :  scale factors state and control
%                 la    :  co-state
%                 maf   :  max(abs(f(x,u)))
%
%           Output:
%                 LQG problem data: A,B,Q,R,V
%
% Programmer: GvW26-1-2003
  
  dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;
  if nargout~=6 || nargin<6; error('LQGABQRV requires 5 outputs and at least 6 inputs'); end;

% Determine state and state and control deviations for the call to linvfxu
% Note that the external input perturbations are read from the final 
% dimd positions of the state vector

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Extended state and control indices
%  xind=1:dimx; uind=dimxd+2:dimxd+dimu+1;
  xind=1:dimxd; uind=dimxd+2:dimxd+dimu+1;

% Extended initial state
%  x=[x; zeros(dimd,1)];
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
% State and control deviations for linxfxu
  xdev=xudev(xind); udev=xudev(uind);

% Compute the linearized dynamics about the optimal trajectory excluding the running costs  
  [A,B]=linvfxu(sysdyn,t,x,u,[],xdev,udev,1:dimxd,xind,varargin{:});
  
% Determine Q(t),R(t),M(t) based on the second order terms of f
  Q=zeros(dimxd); V=zeros(dimxd); R=zeros(dimu); M=zeros(dimxd,dimu);
  [qmr]=slqmr(sysdyn,dims,t,x,u,xudev,ds,'f',la,varargin{:});

% frq > 1: Feedback gain increases
% frq < 1: Feedback gain reduces
%  frq=1;
%  qrf=diag([1/frq*ones(1,dimxd+1) ones(1,dimu)]);
%  qmr=qrf*qmr*qrf;

% Obtain Q,M,R from partitioned matrix
  Q(xind,xind)=qmr(xind,xind); R=qmr(uind,uind);
  M(xind,1:dimu)=qmr(xind,uind);
  
% Obtain V and the matrices to be added to Q,R 
  [V,dum,dum,dum,Qadd,Radd]=feval(lqgmatr,dimxd,dimu,1,1);
  Q=Q+Qadd; R=R+Radd;
