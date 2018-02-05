% readme.m    : Digital optimal control & LQG compensation including
%               a temporal and differential stabilizability & detectability
%               analysis. Presented in the paper:
%               L.G. Van Willigenburg, W.L. De Koning, 2013,
%               "Temporal and one-step stabilisability and detectability of
%                discrete time linear systems",
%               IET Control Theory & Applications, 7, 1, 151-159
%               Download this paper from:
%               http://gvw007.yolasite.com/resources/TempOneIET2013.pdf
%
% Main m-files
% exall.m     : Script that runs ex1 and ex2 the two paper examples
% tmpstdt.m   : Function performs the temporal and differential 
%               stabilizability and detectability analysis.
% 
% Additional problem specific m-files
% ex1.m       : Script example 1. Script solves a continuous optimal
%               control problem with a bang-bang control and then
%               uses the solution to initialize the associated digital 
%               optimal control problem that is solved by opconstr1.
%               Next the digital LQG compensator is computed and
%               the temporal and differential stabilizability and
%               detectability analysis is performed.
% tlsbang.m   : Script uses constr to compute the open loop bang-bang optimal control
% tlsdop.m    : Function that specifies the digital optimal control problem data
% tlsdyn.m    : function that contains the dynamics in the Mayer formulation.
%               format [dxdt]=f(t,x,u,varargin)
% tlsout.m    : function that contains the system output
%               format [y]=g(t,x,u,varargin)
% tlsphi.m    : function that contains the terminal costs
%               format [phi]=phi(t,xf)
%
% ex2.m       : Main script example 2; Solves a digital optimal control
%               problem by opconstr1. Next the digital LQG compensator is
%               computed and the temporal and differential stabilizability
%               and detectability analysis is performed.
% tlsdop1.m   : Function that specifies the digital optimal control problem data
% tlsdyn1.m   : function that contains the dynamics in the Mayer formulation.
%               format [dxdt]=f(t,x,u,varargin)
% tlsout1.m   : function that contains the system output
%               format [y]=g(t,x,u,varargin)
% tlsphi1.m   : function that contains the terminal costs
%               format [phi]=phi(t,xf)
%
% Additional general m-files:
% opconstr1.m : Function computes the digital open loop optimal control
% lqgcomp.m   : Function computes the digital LQG compensator
% tmpstdt.m   : Function performs the temporal and differential stabilizability and
%               detectability analysis.
% stair.m     : Draws staircase functions that represent digital controls
%
% Gerard van Willigenburg / Willem de Koning 17-3-2015
help readme;