function [sx]=signsc(x,left,right,beta)
% SIGNSC: Scaled sign function
% 
%         [sx]=signsc(x,left,right,beta)
%
%         signsc(x)=left if sign(x)=-1
%         signsc(x)=right if sign(x)=1
%         signsc(x)=(left+right)/2 if sign(x)=0
%         beta : small positive scalar used for smoothing (optional)
%                Smoothing is only applied if beta is specified
%
%         x,left,right must have identical dimensions
%         They are treated elements wise.
%
% GvW 21-1-2010

% Check inputs
if nargin<3 || nargin>4; error(' 3 or 4 inputs required'); end
[n,m]=size(x);
[n1,m1]=size(right);
if n1~=n || m1~=m;
  error(' right must have the same size as x');
end
[n1,m1]=size(left);
if n1~=n || m1~=m;
  error(' left must have the same size as x');
end
if nargin==4;
  [n1,m1]=size(beta);
  if n1~=1 || m1~=1 || beta<=0
    error(' beta must be a positive small scalar');
  end
end

% Computation
if isa(x,'gradient') || isa(x,'hessian')
  sx=(x>0).*right+(x<0).*left+0.5*(x==0).*(left+right);
elseif nargin==4
  av=0.5*(right+left); amp=0.5*(right-left);
  sx=av+amp.*x./sqrt(x.*x+beta*beta);
else
  sx=right+0.5*(sign(x)-1).*(right-left);
end