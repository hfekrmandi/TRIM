function [x0m,G,H]=lqgx0gh(phi,lqgmatr,dims,tf,xf,xudev,ds)
% LQGX0GH:   LQG compensation problem specification: x0,G,H
%            
%            [x0m,G,H]=lqgx0gh(phi,lqgmatr,dims,tf,xf,xudev,ds)
%
%            Input:
%              dims : dimensions
%           lqgmatr : Name of file containing the LQG problem data
%              tf   : final time
%              xf   : Final state excluding running costs
%              ds   : Scaling factors
%
%            Output:
%              x0m  : mean x0 excluding running costs
%              G,H  : Initial state intensity and final state penalty matrix
%
% Programmer: GvW 26-1-2003

dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;

if nargout~=3 || nargin~=7; error('LQGX0GH requires 6 input and 3 outputs'); end;

[wp,meth,mt]=lqgwfgi('phi');

% State indices
xind=1:dimxd;

[H]=slqmr(phi,dims,tf,xf,[],xudev(xind),ds,'phi',[],[]);
H=wp*H; x0m=zeros(dimxd,1);

% Obtain x0m, G and the matrix to be added to H
[V,W,G,x0m,Qadd,Radd,Hadd]=feval(lqgmatr,dimxd,dimu,1,1); H=H+Hadd;

% % Extended state indices
% xind=[1:dimxd+1];
% 
% % Extended initial state covariance
% dsh=ds(xind); G=diag(1e-6*ones(1,dimxd)./dsh'./dsh');