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

function [dx]=tlsdyn1(t,x,u,varargin)

% Parameters
aalpha = 0.01227; bbeta = 0.145e-3; c = 2060; g0 = 9.81;

% state and control dependent expressions
D = aalpha*x(2)*x(2)*exp(-bbeta*x(1)); g = g0;

% state derivatives
dx = [x(2); (u(1)*c-D)/x(3)-g; -u(1); 0];