% TLSDYN: Dynamics system with bang-bang optimal control
%
%         Input:
%         t      : current time
%         x      : current state
%         u      : current input
%         p      : parameter vector
%         td     : external inputs [t d], t: columnvector d: stacked row vectors
%
%         Output:
%         dx     : dx/dt
%
% GvW 28-8-2007

function [dx]=tlsdyn(t,x,u,varargin)
p=0; dx=[-x(1)-x(2)+u(1)+0.05*u(1)*x(1); x(1)*(x(2)-p(1)); (x(2)-p(1))*(x(2)-p(1))];