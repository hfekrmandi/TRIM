function [xe,ye]=plotellip(A,b,c,linetype)
% [xe,ye]=plotellip(A,xc,linetype);
%   plots ellipsoid (x-xc)*A*(x-xc) = 1 in R^2, or
% [xe,ye]=plotellip(A,b,c,linetype);
%   plots ellipsoid x'*A*x + b'*x + c = 0 in R^2

if nargin==2
  type1=1;  xc=b;  linetype='-';
elseif nargin==3 & isstr(c)
  type1=1;  xc=b;  linetype=c;
elseif nargin==3
  type1=0;  linetype='-';
else
  type1=0;
end

A=.5*(A+A');
if min(eig(A))<=0
  disp('invalid ellipsoid -- A not positive definite'), xe=[];, ye=[];, return
end
if ~type1
  xc=A\(-b./2);
  k=xc'*A*xc-c;
  if k<=0
    disp('invalid ellipsoid'), xe=[];, ye=[];, return
  end
  A=A/k;
end

x0 = xc(1); y0 = xc(2);
R = chol(A);
theta = linspace(-pi,pi,100);
xy_tilde = [cos(theta); sin(theta)];
invR = inv(R);
xy_bar = invR*xy_tilde;
xe = xy_bar(1,:)+x0;
ye = xy_bar(2,:)+y0;
plot(xe,ye,linetype);

