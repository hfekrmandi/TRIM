function [f1,f2,tt,xt,yt]=tlsdop1(u,x,t,ts,flg,varargin)
% TLSDOP  : Example of a digital optimal control problem with fixed final time
%           and terminal equality/inequality constraints PSI and bounded inputs
%           in the Mayer formulation i.e. L(x,u,t) is last component of f(x,u,t)
%           and so the last component of x is integral L(x,u,t)dt. Also J=PHI.
%
%           [f1,f2]=tlsdop(u,x,t,ts,flg,varargin)
%
% Inputs :
%           u	       : constant control over current sampling interval (mx1)
%           x	       : initial state (flag=1) final state (flag=2)
%           t	       : current sampling time
%           ts	     : length current sampling interval
%           flg      : flag
%        varargin    : additional inputs: varargin{1:4}={[tu] [x0] [pv] [extinp]}
%           tu       : initial guess solution: [t u] table
%           pv       : parameter vector system
%           extinp   : external input trajectories: [t u] table
%
% Outputs :
%          flg=1    : f1=final state
%          flg=2    : f1=[PHI;PSI], f2=d[PHI;PSI]/dx, J=PHI (because of the Mayer formulation)
%          flg=3    : f1=df/dx, f2=df/du
%          flg=4    : f1=Sampling instants final time and initial control matrix
%                     f2=fixed stepsize numerical integration
%                     f1=[[t0;t1;..;tN],[u0;u1;..uN]],f2=stepsize numerical integration
%                     Note that uN does not influence the performance but must be specified
%          flg=5    : f1=Initial state
%                     f2=System parameter vector
%          flg=6    : f1=Control bounds [ulow,uhigh]
%                     f2=Maximum absolute control change during a single iteration
%                     Either f1 or f2 must be specified (non-empty). If only f2
%                     is specified f1 is assumed +/- 10 times f2.
%          flg=7    : f1=State inequality constraints PSIT and associated derivatives f2=dPSIT
%                     Specify f1=[], f2=[] if not present
%          flg=8    : f1=Number of equality constraints that appear first in PSI
%                     Specify f1=0 if not present
%          flg=9    : f1=0 or [] then dHdu used for gradient computation, f1~=0 then dHdu not used
%
% GvW 11-02-2010

dt=0.1; % Fixed stepsize numerical integration
f1=[]; f2=[]; tt=[]; xt=[]; yt=[]; % Output empty matrices unless specified
if flg==1
  % Integrate system equations from x(t) to x(t+ts), f1=x(t+ts)
  [tt,xt,yt]=ode23fs('tlsdyn1','tlsout1',[t t+ts],x,dt,u,varargin{3:end});
  f1=xt(end,:)';
elseif flg==2
  % Specify the terminal objective J=PHI and the terminal equality/inequality
  % constraints PSI=0, PSI<0: f1=[PHI;PSI], f2=d[PHI;PSI]/dx
  % Equality constraints must appear first in PSI, see also flg==8!
    f1=[-x(1); -x(3)+67.9833]; f2=[-1 0 0 0; 0 0 -1 0];
%     phipsi=@(x)([tlsphi1(0,x); -x(3)+67.9833]);
%     f12=phipsi(gradientinit(x)); f1=f12.x; f2=f12.dx;
elseif flg==3
% Specify the Jacobian df/dx (includes dL/dx) f1 and df/du (includes dL/du) f2
  [f1,f2]=dlfxuad('tlsdyn1',t,ts,x,dt,u,varargin{3:end});
elseif flg==4
  % Specify matrix with sampling instants final time and initial control f1
  % and fixed time step numerical integration f2
  % f1 is overruled if an input file is specified in the call to opconstr
  if ~isempty(varargin); tu=varargin{1}; f1=tu; % External initial guess
  else
    tf=56; N=28; Ts=tf/N; tsar=[0:Ts:N*Ts]'; u=zeros(size(tsar));
    u(tsar<=4)=9.525515; %u(tsar>4 & tsar<=30)=2.13; 
    th=tsar(tsar>4 & tsar<=46); th=th-th(1);
    u(tsar>4 & tsar<=46)=2.13+(4.19-2.13)*(th/th(end)).^6; u(tsar>46)=0;
    %u(tsar>4 & tsar<=46)=4.19; u(tsar>46)=0;
    f1=[tsar u]; stair(tsar,u); title('Initial control');
    xlabel('t'); ylabel('u'); drawnow;
  end
  f2=dt;
elseif flg==5
  % Specify initial state f1 and and the system parameter vector f2
  if length(varargin)>1; x0=varargin{2};
  else x0=[0;0;214.839;0]; end; f1=x0;
  if length(varargin)>2; pv=varargin{3}; else pv=[]; end; f2=pv;
elseif flg==6
  % Specify matrix with control bounds f1 and vector with maximum control
  UL=[0]; UU=[9.525515]; u_constr=[UL UU]; f1=u_constr; f2=0.5*(UU-UL);
elseif flg==7
  % Specify state inequality constraints PSIT, dPSIT
  f1=[-x(1); -x(2); -x(3)+67.9833; x(3)-214.839];
  f2=[-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 1 0];
elseif flg==8
  % Specify number of equality constraints that appear FIRST in PSI
  f1=0;
elseif flg==9;
  f1=0;
end