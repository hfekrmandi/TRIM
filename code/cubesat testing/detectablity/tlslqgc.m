function [V,W,G,x0m,Qa,Ra,Ha,nc,Wf,Wg,WI,meth,mt,nmc]=tlslqgc(n,m,l,ns)
% TLSLQGC  : LQG compensator design matrices Q,R,H,V,W,G,x0m and state dimension nc
%            
%
%           Input:
%                 n  :  dimension state
%                 m  :  dimension control
%                 l  :  dimension output
%                 ns :  number sampling instants
%
%           Output:
%                 V,W,G,x0m      : LQG problem data
%                 Qa,Ra,Ha : LQG problem data added to Q, R obtained from
%                                  the second-order non-linear terms in f and the second-order
%                                  terms in the costs J.
%                 nc(i), i=1,ns  : Prescibed dimensions of the LQG compensator state
%                                  at each sampling instant i=1,..,ns. If a
%                                  single value is specified this value
%                                  will be used for each sampling instant.
%                 Wf             : Wf, possibly a vector with different values for Wf
%                                  the weighting factor of the second-order terms in f
%                 Wg             : Wg, possibly a vector with different values for Wg
%                                  the weighting factor of the second-order terms in g
%                 WI             : WI, possibly a vector with different values for WI
%                                  the weighting factor of the second-order terms in J
%                 meth           : 0 then Wf, WI redundancy incorporated
%                 mt             : 1,2,3 method of approximating indefinite quadratic forms see cin2nn
%                 nmc            : Number of runs each Monte Carlo simulation
%                                
% Programmer: GvW 26-1-2003

% LQG design matrices
V=eye(n); W=eye(l); G=eye(n); x0m=zeros(n,1);

% Matrices to be added
Qa=eye(n); Ra=eye(m); Ha=eye(n);

% Compensator state dimension
nc=n; % Full-order LQG compensator

% Weighting factor vectors for wf, wg, wp
Wf=0; Wg=0; WI=0;
%WI=[1 2 5 7.5 10 12.5];

% Computation methods (redundancy included, semi-positive diagonal matrices)
meth=0; mt=2;

% Number of runs 1 Monte Carlo simulation
nmc=1;