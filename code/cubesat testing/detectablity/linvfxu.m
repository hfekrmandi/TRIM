function [A,B]=linvfxu(fcn,t,x,u,para,xpert,upert,ind,xind,varargin)
% LINVFXU : Differentiate a vector function fcn of x,u
%           to x, result A, and u, result B
%
%           [A,B]=linvfxu(fcn,t,x,u,para,xpert,upert,ind,xind,varargin)
%
%           para : optional parameter to calculate perturbations
%           xpert=para(1)+1e-3*para(1)*abs(x);
%           upert=para(1)+1e-3*para(1)*abs(u);
%           para(1) : default 1e-5;
%           xpert : perturbations x for finite difference computation
%                   The dimension n of xpert may be smaller than n of x
%                   In that case fcn is differentiated only w.r.t.
%                   the first n components of x. If xind is specified
%                   xpert are perturbations associated to x(xind)
%           upert : perturbations u for finite difference computation
%           ind   : optional array of indexes for the function fcn
%                   If specified only the corresponding components of fcn
%                   are differentiated respectively
%           xind  : optional array of indexes for the perturbations xpert
%                   If specified only the corresponding components of x
%                   are varied in the specified order respectively
%
%           Format f(x,u):
%           fvalue = f(t, x, flag, u);
%           in which flag is not used by linvfxu

%
if ~(ischar(fcn) || isa(fcn,'function_handle')); 
  error(' fcn must be a string or function handle'); 
end

[nx,nh]=size(t);
if nx~=1 || nh~=1; error('  t must be a scalar'); end;

[nx1,nh]=size(x);
if nh~=1; error('  x must be a column vector'); end;

[nu,nh]=size(u); if max(nu,nh)==0; uy=0; else uy=1; end;
if uy && nh~=1; error('  u must be a column vector'); end;

if nargin<8;  ind=1:nx1; end;
if isempty(ind); ind=1:nx1; end
  
if nargin<9;  xind=1:min(nx1,max(ind)); end;
if isempty(xind); xind=1:min(nx1,max(ind)); end

if nargin<5; para(1)=1e-5; end;
if isempty(para); para(1)=1e-5; end;

if nargin<6
  xpert=para(1)+1e-3*para(1)*abs(x(xind));
end

if uy && nargin<7
  upert=para(1)+1e-3*para(1)*abs(u);
end

if isempty(xpert)
  xpert=para(1)+1e-3*para(1)*abs(x(xind));
end

if uy && isempty(upert)
  upert=para(1)+1e-3*para(1)*abs(u);
end

[nx,nh]=size(upert);
if uy && nh~=1; error('  upert must be a column vector'); end;

[nx,nh]=size(xpert);
if nh~=1; error('  xpert must be a column vector'); end;
if nx>nx1; error(' dimension xpert must be <= dimension x'); end;

in=0;
if nargin>=8 ;
  if ~isempty(ind); in=1;
  else
    [nh,mh]=size(ind); if nh~=1 || mh~=1; error('  ind must be a vector'); end;
  end
end

xin=0;
if nargin>=9 ;
  if ~isempty(xind); [nh,mh]=size(xind);
    if nh~=1 && mh~=1; error('  xind must be a vector'); end;
    xin=1;
  end
end
ix=1:nx1; iu=nx1+1:nx1+nu;

% Gradients (intlab)
if exist('gradientinit','file')
  xu=gradientinit([x;u]);
  dx = feval(fcn, t, xu(ix), xu(iu), varargin{:}); AB=full(dx.dx);
  A=AB(:,ix); B=AB(:,iu);

  % Adjustment
  if in && xin
    A=A(ind,xind); B=B(ind,:);
  elseif in
    A=A(ind,:); B=B(ind,:);
  elseif xin
    A=A(:,xind);
  end

elseif exist('fmad','file')% Gradients (tomlab)
  xu=fmad([x;u],eye(nx1+nu));
  dx = feval(fcn, t, xu(ix), xu(iu), varargin{:});
  AB=full(getinternalderivs(dx));
  A=AB(:,ix); B=AB(:,iu);

  % Adjustment
  if in && xin
    A=A(ind,xind); B=B(ind,:);
  elseif in
    A=A(ind,:); B=B(ind,:);
  elseif xin
    A=A(:,xind);
  end
  
else
  %  Gradients finite differences
  dx = feval(fcn, t, x, u, varargin{:});
   
  if in; 
    if any(rem(ind,1))~=0 || min(ind)<1
      error(' ind must contain positive integers');
    end;
    if max(ind)>size(dx,1); error(' max(ind) > size(fcn)'); end;
    dx=dx(ind);
  end;
  olddx=dx; oldx=x; oldu=u;

  if xin; 
    if size(xpert,1)~=length(xind);
      error(' xind and xpert must have the same length');
    end;
    if any(rem(xind,1))~=0 || min(xind)<1
      error(' xind must contain positive integers');
    end;
    if max(xind)>size(x,1); error(' max(xind) > size(x)'); end;
  else
    xind=1:size(xpert,1);
  end

  nx=size(xpert,1); nf=max(size(dx)); A=zeros(nf,nx); B=zeros(nf,nu);

  % A matrix
  j=1;
  for i=xind;
    x(i)=x(i)+xpert(j);
    dx = feval(fcn, t, x, u, varargin{:});
    if in; dx=dx(ind); end;
    A(:,j)=(dx-olddx)./xpert(j); j=j+1;
    x=oldx;
  end

  if uy && nargout>1
    % B matrix
    for i=1:nu;
      u(i)=u(i)+upert(i);
      dx = feval(fcn, t, x, u, varargin{:});
      if in; dx=dx(ind); end;
      B(:,i)=(dx-olddx)./upert(i);
      u=oldu;
    end
  end
end
