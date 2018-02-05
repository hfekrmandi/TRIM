function [C,W]=lqgcw(sysout,lqgmatr,dims,t,x,u,xudev,ds,l,varargin)
% LQGCW:   LQG compensator design matrices C,W
%            
%          [C,W]=lqgcw(sysout,lqgmatr,dims,t,x,u,xudev,ds,l,varargin)
%
%          Input:
%               sysout : String with name of function containing the system output
%             lqgmatr  :  Name of file containing the LQG problem data
%                 dims :  dimensions
%                 t    :  current time
%                 x    :  current state including running costs (Mayer formulation)
%                 u    :  current control
%                 ds   :  scale factors state and control
%                 dsd  :  scale factors external inputs
%                 l    :  dimension output
%            varargin  :  additional inputs. varargin{1:2}={[pv] [extinp]}
%
%          Output:
%                 LQG problem data: C(tk), W(tk)
%
% Programmer: GvW 26-1-2003

  dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;
 
  if nargout~=2 || nargin<9; error('LQGCW requires at least 9 inputs and 2 outputs'); end;
  
% Determine state and state and control deviations for the call to linvfxu
% Note that the external input perturbations are read from the final 
% dimd positions of the state vector

% Extended state control and output indices
  xind=1:dimxd; uind=dimxd+2:dimxd+dimu+1; yind=dimxd+dimu+2:dimxd+dimu+dimy+1;

% Extended initial state
  %x=[x; zeros(dimd,1)];
  
% State and control deviations for linvxfxu
  xdev=xudev(xind); udev=xudev(uind);

% Compute the linearized output about the optimal trajectory excluding the running costs
  [C,dh]=linvfxu(sysout,t,x,u,[],xdev,udev,1:dimy,xind,varargin{:});
  
% Obtain W matrix
  [dum,W]=feval(lqgmatr,dimxd,dimu,dimy,1);
  
% % Specify the output noise (uncertainty)
%   W=diag(1e-6*ones(dimy,1)./ds(yind)./ds(yind));