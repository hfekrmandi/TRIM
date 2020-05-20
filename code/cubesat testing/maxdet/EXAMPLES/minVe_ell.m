function [A,b,c] = minVe_ell(As,bs,cs)
% [A,b,c] = minVe_ell(As,bs,cs);
%
% find the minimum-volume ellipsoid in R^2
%   { x | x^T*A*x + 2*b^T*x + c < 0 }
% containing given ellipsoids
%   { x | x^T*A_i*x + 2*b_i^T*x + c_i < 0 }, i = 1,...,L
% 
% As = [A1 A2 ... AL]:  L 2-by-2 positive definite matrices
% bs = [b1 b2 ... bL]:  L 2-vectors
% cs = [c1 c2 ... cL]:  L scalars
%
% NOTE: minVe_ell gives a demo if no input argument is given.

% last major update: 01/11/96

if nargin<1
As=[0.1355  0.1148  0.6064 -0.1022  0.7127 -0.0559  0.2706 -0.1379  0.4008 -0.1112
    0.1148  0.4398 -0.1022  0.7344 -0.0559  0.9253 -0.1379  0.2515 -0.1112  0.2107];
bs=[-0.2042    0.8259   -0.0256    0.1827    0.3823
     0.0264   -2.1188    1.0591   -0.3844   -0.8253];
cs=[ 0.2351    5.8250    0.9968   -0.2981    2.6735];
end

L = size(bs,2);    % L ellipsoids given

% check ellipsoids
for i=1:L
  A=As(:,2*i-1:2*i); b=bs(:,i); c=cs(i);
  eigA=eig(A);
  if (~isreal(eigA) | min(eigA)<=0 | c>=b'*inv(A)*b)
    disp([int2str(i) 'th ellipsoid is invalid'])
    A=[]; b=[]; c=[]; return
  end
end

% construct F
%
% diag( [  A, ... 
%          tau_i[Ai bi 0; bi' ci 0; 0 0 0] - [A b 0; b' -1 b'; 0 b -A],
%          tau_1, ..., tau_L )

% number of variables: 
novars = L + 5;   % A, b, tau1 ... tauL
 
% block 1...L
for i = 1:L
F = [F; ...
      [ 0;  0;  0;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  1;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  0;  0;  0], ...   % const
      [-1;  0;  0;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  0;  1;  0;
        0;  0;  0;  0;  0], ...   % A1
      [ 0; -1;  0;  0;  0;
       -1;  0;  0;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  0;  0;  1;
        0;  0;  0;  1;  0], ...   % A2
      [ 0;  0;  0;  0;  0;
        0; -1;  0;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  0;  0;  0;
        0;  0;  0;  0;  1], ...   % A3
      [ 0;  0; -1;  0;  0;
        0;  0;  0;  0;  0;
       -1;  0;  0; -1;  0;
        0;  0; -1;  0;  0;
        0;  0;  0;  0;  0], ...   % b1
      [ 0;  0;  0;  0;  0;
        0;  0; -1;  0;  0;
        0; -1;  0;  0; -1;
        0;  0;  0;  0;  0;
        0;  0; -1;  0;  0], ...   % b2
       zeros(25,i-1), ...
      reshape([As(:,2*(i-1)+1:2*(i-1)+2), bs(:,i), zeros(2,2);
               bs(:,i)', cs(i), 0, 0;
               zeros(2,5)],25,1), ...    %tau
       zeros(25,L-i)];
end;

% last block
F = [F; zeros(L,6), eye(L,L)];
F_blkszs = [5*ones(1,L) ones(1,L)];

% forming G
G = [ [0; 0; 0; 0], ...   % const
      [1; 0; 0; 0], ...   % A1
      [0; 1; 1; 0], ...   % A2
      [0; 0; 0; 1], ...   % A3
      zeros(4,L+2)];      % b1, b2 and taus
G_blkszs = 2;

% phase1
[x0,Z,W,ul,infostr]=phase1(F,F_blkszs,G,G_blkszs,100,1e-3,1e-3,100);

% check feasibility
if ~strcmp('feasible',infostr)
  disp(infostr)
  disp('problem is not feasible!')
  A=[];, b=[];, c=[];, return
end

% call maxdet
c = zeros(novars,1);
[x,Z,W,ul,hist,infostr]=maxdet(F,F_blkszs,G,G_blkszs,c,x0,...
                               zeros(size(F,1),1),zeros(size(G,1),1),...
                               1e-5,1e-5,100,100);

% form the ellipsoid
A=[x(1) x(2); x(2) x(3)];
b=[x(4); x(5)];
c=-1+b'*inv(A)*b;

% make plot
clg, hold on
for i=1:L
  plotellip(As(:,2*(i-1)+1:2*(i-1)+2),2*bs(:,i),cs(i));
end
[tmpx,tmpy]=plotellip(A,2*b,c,':');
tmp=ceil(1.25*max(max(abs([tmpx,tmpy]))));
axis('square')
axis('off')
hold off