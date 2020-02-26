function edorp(sysdyn,sysout,phi,lqgmatr,uystar,edorp)
% EDORP :  Function that computes the Equivalent discrete optimal LQG regulator problem
%
%          edorp(sysdyn,sysout,phi,lqgmatr,uystar,edorp)
%
%          Input:
%          sysdyn : String with name of function containing the system dynamics
%                   format function [dxdt]=f(t,x,u,varargin)
%          sysout : String with name of function containing the system dynamics
%                   format function [y]=g(t,x,u,varargin)
%          phi    : String containing the name of phi the terminal costs
%                   format function [phi]=phi(t,x(tf))
%        lqgmatr :  Name of file containing the LQG problem data
%          uystar : String with name of file containing the optimal control data
%
%          Output:
%          edorp  : String containing the name of the file that contains the output data
%
%          See also lqgedat.m intmdyn.m
%
%          L.G. van Willigenburg 18-6-'03.
%
%  clear; close all;
  
% Check inputs
  if nargin~=6; error(' 6 input required'); end
  if ~isstr(sysdyn); error(' 1st input must be a string'); end
  if ~isstr(sysout); error(' 2nd input must be a string'); end
  if ~isstr(phi); error(' 3rd input must be a string'); end
  if ~isstr(lqgmatr); error(' 4th input must be a string'); end
  if ~isstr(uystar); error(' 5th input must be a string'); end
  if ~isstr(edorp); error(' 6th input must be a string'); end

% Load optimal control data
  load(uystar)

% Compute Equivalent discrete optimal LQG regulator problem
  satimes=tk'; sa1=satimes(1:end-1); sa2=satimes(2:end);
  nit=[round((sa2-sa1)/dt) 1]; satimes=[satimes; nit];
  optcon=[tk ustar xstar lk];
  
  [lqgarr,costed,tt,qt,rt,mct,lat,vt,wt]=lqgedat('lqgdsamp',sysdyn,sysout,phi,lqgmatr,dims,satimes,optcon,xudev,ds,maf,pv,extinp);

% Save equivalent discrete LQG problem
  save(edorp); disp('EDORP saved');