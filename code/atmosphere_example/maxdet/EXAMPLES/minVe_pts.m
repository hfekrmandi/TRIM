function [A,b,c] = minVe_pts(X)
% [A,b,c] = minVe_pts(X);
%
% find the minimum-volume ellipsoid in R^2
%      { x | x^T*A*x + 2*b^T*x + c < 0 }
% containing given points xi, i=1,...,n
%
% X:  2-by-n matrix contains n given points in R^2. 
%     X = [x1 ... xn]
%
% NOTE: minVe_pts gives a demo (requires exp_data.mat) if no input 
%       argument is given

if nargin<1
  load exp_data
  X=V;
end
n=size(X,2);

% form G
G=[0 0 0 0;
   1 0 0 0;
   0 1 1 0;
   0 0 0 1;
   0 0 0 0;
   0 0 0 0]';
G_blkszs=2;

% form F
F=[];  F_blkszs=[];
for i=1:n
  F=[F; [1 0 0 0 1 0 0 0 1;
         0 0 X(1,i) 0 0 0 X(1,i) 0 0;
         0 0 X(2,i) 0 0 X(1,i) X(2,i) X(1,i) 0;
         0 0 0 0 0 X(2,i) 0 X(2,i) 0;
         0 0 1 0 0 0 1 0 0;
         0 0 0 0 0 1 0 1 0]'];
  F_blkszs=[F_blkszs 3];
end

% find feasible initial point
radinv = min(.5/sqrt(max(X(1,:).^2+X(2,:).^2)),.01);
x0=[radinv; 0; radinv; 0; 0];

% call maxdet
[x,Z,W,ul,hist,infostr]=maxdet(F,F_blkszs,G,G_blkszs,zeros(5,1),x0,...
           zeros(size(F,1),1),zeros(size(G,1),1),1e-5,1e-5,100,100);

% form the ellipsoid
E=[x(1) x(2); x(2) x(3)];
d=[x(4); x(5)];

A=E'*E;
b=E'*d;
c=d'*d-1;

% make plot
clg
plot(X(1,:),X(2,:),'+')
hold on
plotellip(A,2*b,c,':');
tmp=ceil(1.25*max(max(abs(X))));
axis([-tmp tmp -tmp tmp])
axis('square')
axis('off')
hold off
