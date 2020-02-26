function [nc]=lqgnc(lqgmatr,dims,ns,x0)
% LQGNC:   LQG compensator dimension specification
%            
%          [nc]=lqgnc(lqgmatr,dims,ns,x0)
%
%          Input:
%         lqgmatr : name fof file containing the LQG problem data
%            dims : dimensions
%            ns   : number of sampling intervals
%            x0   : initial state
%
%          Output:
%            nc  : dimensions of the compensator state
%
% Programmer: GvW 26-1-2003

dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;
 
if nargin~=4; error(' LQGNC requires 4 input arguments'); end;
if nargout~=1 error('LQGNC requires 1 output argument'); end;

% Get compensator state-dimension
[dum,dum,dum,dum,dum,dum,dum,nc]=feval(lqgmatr,dimxd,dimu,dimy,ns);

% Compute max-min dimensions
if nc==dimxd || ~any(nc-dimxd);
  nc=nc*ones(1,ns+1);
else
  if max(abs(x0))==0; nc0=0; else nc0=1; end;
  h1=nc0:dimu:dimu*ns+nc0; h2=fliplr(0:dimy:dimy*ns); h3=nc*ones(1,ns+1);
  nc=min(h1,h2); nc=min(nc,h3);
end