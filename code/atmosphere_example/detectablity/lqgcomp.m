function lqgcomp(sysdyn,sysout,phi,lqgmatr,infile,mro,opt)
% LQGCOMP: LQG compensator computation and Monte Carlo simulation
%          The compensator is saved in file comp.mat.
%          The Monte Carlo simulation results are saved in mocarlo.mat.
% 
%         lqgcomp(sysdyn,sysout,phi,lqgmatr,infile,mro,opt)
%
%         Input:
%         sysdyn  : String with name of function containing the system dynamics
%                   format function [dxdt]=f(t,x,u,varargin)
%         sysout  : String with name of function containing the system dynamics
%                   format function [y]=g(t,x,u,varargin)
%         phi     : String containing the name of phi the terminal costs
%                   format function [phi]=phi(t,x(tf))
%         lqgmatr : String with name of file containing the LQG problem data
%         infile  : String with the file name of the file containing
%                   the solution of the open-loop digital optimal control
%                   This file can be obtained from opconstr 
%         mro     : Method to compute reduced-order LQG compensator
%                   1 : Optimal Projection equations (default), function dprotva
%                   2 : Lyapunov equations, function dprotvla
%                   3 : UDU factored Lyapunov equations, function dprtvlau
%         opt     : options for
%                   1 : dprotva
%                   2 : dprotvla
%                   3 : dprtvlau
%
%         Output files:
%         comp.mat    : digital optimal compensator Matlab 4 format
%         mocarlo.mat : Monte Carlo simulation results
%         edorp.mat   : Equivalent discrete-time LQG problem data file
%         uystar.mat  : digital optimal control and output
%         uy.dat      : digital optimal control and output ascii format
%         digcomp.mat : digital optimal compensator and associated data
%         xfkl.dat    : digital optimal compensator ascii format
%
%         Example: lqgcomp('helidyn','heliout','heliphi','helilqgc','dopco')
%
% GvW 24-01-2010

% Check input
if nargin<5; error(' At least 5 inputs required'); end
if ~(ischar(sysdyn) || isa(sysdyn,'function_handle')); 
  error(' 1st argument must be a string or function handle'); 
end
if ~(ischar(sysout) || isa(sysout,'function_handle')); 
  error(' 2nd argument must be a string or function handle'); 
end
if ~(ischar(phi) || isa(phi,'function_handle')); 
  error(' 3rd argument must be a string or function handle'); 
end
if ~ischar(lqgmatr); error(' 4th input must be a string'); end
if ~ischar(infile); error(' 5th input must be a string'); end
if ~exist('mro','var'); mro=1; end;
if ~exist('opt','var'); opt=[]; end;

global W_f W_g W_I me_th m_t;

% Get the different values of wf, wg, wp to be evaluated through Monte Carlo simulation
% and nmc the number of runs for each Monte Carlo simulation
[dum,dum,dum,dum,dum,dum,dum,dum,Wfv,Wgv,WIv,me_th,m_t,nmc]=feval(lqgmatr,1,1,1,1);

% Get computational method
[WfgI,meth,mt]=lqgwfgi('f');

% Initialization
WfgIt=[]; fail=[]; failt=[]; acdev=[]; acdevt=[]; critt=[]; xmx=[]; loop=0;

format short e; format compact;

for W_f=Wfv
  for W_I=WIv
    for W_g=Wgv

% Simulate the optimal open loop system without disturbances
% for verification and to compute the optimal output y(t(k))
loop=loop+1; simopt(sysdyn,sysout,phi,infile,'uystar','uy.dat'); %close all;

% Compute the optimal compensator and store
edorp(sysdyn,sysout,phi,lqgmatr,'uystar','edorp');
diglqg(lqgmatr,'edorp','digcomp','xfkl.dat','comp',opt); 
%close all;
mulnoi=1; costrt=[]; kk=0;

% Monte Carlo simulations of the associated closed loop systems
flag=0; save 'mcsim.txt' flag -ascii;
hw=waitbar(0,['Monte Carlo simulation : ' num2str(nmc) ' runs']);
for sd=1:nmc
  kk=kk+1; costrt(kk,1)=sd;
  rp=0; cw=1; simcsys; %close all;
  costrt(kk,2)=costr; xmx(sd,:)=max(xmk);
  waitbar(sd/nmc,hw);
%  disp(costrt(end,:));
  if fileflag('mcsim.txt');
    close(hw);
    error(' Termination due to mcsim.txt');
  end
end
close(hw);

% Detect, count and exclude Inf and Nan results
% which represent closed loop system failures (due to instability)
point=[]; pointer=1:size(costrt,1);
for kk=size(costrt,1):-1:1
  if max(isinf(costrt(kk,:))) || max(isnan(costrt(kk,:)))
    pointer(kk)=[]; point=[kk point];
  end
end
costrh=costrt(pointer,:);

% Compute and display the results
crit(1:2)=[cost mean(costrh(:,2))];
h=cumsum(isinf(costrt(:,2))+isnan(costrt(:,2)));
fail(1)=h(end);
acdev(1)=mean(abs(costrh(:,2)-cost));

acdevt=[acdevt; acdev]; failt=[failt; fail];
WfgIt=[WfgIt; W_f W_g W_I]; critt=[critt; crit];
pointt{loop}=pointer; costt{loop}=costrt;

format short g;
disp(' '); disp('           Wf           Wg           WI      Delta J     Failures');
disp([WfgIt acdevt failt]);
 
    end
  end
end
eval(['save mocarlo acdevt failt WfgIt critt pointt costt xmx cost loop meth mt sysdyn sysout phi']);