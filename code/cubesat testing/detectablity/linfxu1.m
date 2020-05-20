function [A,B]=linfxu1(fcn,t,x,u,para,xpert,upert,varargin)
% LINFXU : Differentiate a vector function fcn(x,u)
%          to x, result A, and u, result B
%
%          [A,B]=linfxu(fcn,t,x,u,para,xpert,upert,p1,p2,...)
%
%          para : optional parameter to calculate perturbations
%          xpert=para(1)+1e-3*para(1)*abs(x);
%          upert=para(1)+1e-3*para(1)*abs(u);
%          default 1e-3;
%          xpert    : perturbations x
%          upert    : perturbations u
%          p1,p2,.. : additional parameters passed on to fcn

%
if ~isstr(fcn); error('  fcn must be a string with a function name'); end;

[nx,nh]=size(t);
if nx~=1 | nh~=1; error('  t must be a scalar'); end;

[nx,nh]=size(x);
if nh~=1; error('  x must be a column vector'); end;

[nu,nh]=size(u);
if nh~=1; error('  u must be a column vector'); end;

if nargin<5; para(1)=1e-3; end;
if isempty(para); para(1)=1e-3; end;

if nargin<6; xpert=para(1)+1e-3*para(1)*abs(x); end
if isempty(xpert);  xpert=para(1)+1e-3*para(1)*abs(x); end

if nargin<7; upert=para(1)+1e-3*para(1)*abs(u); end
if isempty(upert);  upert=para(1)+1e-3*para(1)*abs(u); end

% Initialization
oldx=x; oldu=u;
dx = feval(fcn, t, x, u, varargin{:});
olddx=dx;

nf=max(size(dx)); A=zeros(nf,nx); B=zeros(nf,nu);

% A matrix
for i=1:nx;
  x(i)=x(i)+xpert(i);
  dx = feval(fcn, t, x, u, varargin{:});
  A(:,i)=(dx-olddx)./xpert(i);
  x=oldx;
end

% B matrix
for i=1:nu;
  u(i)=u(i)+upert(i);
  dx = feval(fcn, t, x, u, varargin{:});
  B(:,i)=(dx-olddx)./upert(i);
  u=oldu;
end
