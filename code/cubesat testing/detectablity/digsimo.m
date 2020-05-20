function [tt,xt,yt]=digsimo(u,name,x0,tsar,varargin)
% DIGSIMO : Digital simulation of system to obtain intersample behavior
%
%          [tt,xt,yt]=digsimo(u,name,x0,tsar,pv,varargin)
% 
%          Input
%          u         : (m x ni) vector of initial controls
%          name	     : name of the model-file (m-file)
%          x0        : (n x 1) state at time t=0
%          tsar      : array with sampling instants and final time
%          varargin  : additional inputs. Varargin{1:2}={[pv] [extinp]}
%
%          Output
%          tt        : time
%          xt        : states
%          yt        : outputs
%
% GvW: 19-04-2007

ns=length(x0); [nc,N]=size(u);
x=zeros(ns,1); x=x0; xt=[]; yt=[]; tt=[];
% Forward sequencing and store x(:,i):
for i=1:N,  
  [x,f2,th,xh,yh]=feval(name,u(:,i),x,tsar(i),tsar(i+1)-tsar(i),1,varargin{:});
  tt=[tt;th]; xt=[xt;xh]; yt=[yt;yh];
end