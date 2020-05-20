  function [xk,yk,tout,xout,yout,maf]=digsim(fxu,dims,tk,uk,x0,dt,gxu,df,varargin)
%
% DIGSIM : Function simulates a continuous time system
%          with a piecewise constant control using ode23
%
%          Input:
%                 fxu   :  string with function name xdot = fxu(t,x,flg,u) V5.1
%                          or xdot = fxu(t,x) V4.2c2
%                 tk    :  column vector with sampling instants
%                 uk    :  piecewise constant control
%                 x0    :  initial state
%                 dt    :  optional prescribed time steps
%                          within sampling interval that determines tout,xout
%                 gxu   :  optional string with output function name
%                          y = gxu(t,x,flg,u)
%                 df    :  If non-empty and non-zero a waitbar is plotted
%                          and termination can be invoked through editing
%                          the file digsim.txt
%           varargin    :  Additional inputs. Varargin{1:2}=[pv extinp]
%
%	   Output:
%                 xk    :  states at the sampling instants
%                 yk    :  outputs at the sampling instants
%                 tout  :  t0:dt:tN
%                 xout  :  states at tout
%                 maf   :  max(abs(f(x,u))
%

  dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5);

% check inputs
  if nargin<5; error('  At least 5 input argumenst required'); end;
  if ~(ischar(fxu) || isa(fxu,'function_handle')); 
    error(' 1st argument must be a string or function handle'); 
  end
  [nt,mt]=size(tk); [nu,mu]=size(uk);
  if mt~=1; error('  2nd input must be a column vector'); end;
  if  nu<nt-1; error('  tk and uk incompatible'); end;
  
  if nargin<6; dt=0; gxu=[]; df=0;
  elseif nargin<7; gxu=[]; df=0;
  elseif nargin<8; df=0; end;
  if ~isempty(gxu)
    flgy=1;
    if ~(ischar(gxu) || isa(gxu,'function_handle')); 
      error(' 5th argument must be a string or function handle'); 
    end
  else flgy=0;
  end
  
  if max(size(dt))~=1; dt=0; end;

%  global u  %V4.2c2

% integrate the system and
% save the states at the sampling instants
  xi=x0; xk=x0'; maf=zeros(size(x0)); tout=[]; xout=[]; yout=[];
  flag=0; save digsim.txt flag -ascii;
  if df; hw=waitbar(0,'Simulation digital control'); end;
  for i=1:nt-1
     ti=tk(i); tf=tk(i+1); th=(tf-ti)/2; u=uk(i,:)';
     if dt<=0; dt=th; else min(dt,th); end; tv=ti:dt:tf;
     [t,x,y,mf]=ode23fs(fxu,gxu,tv,xi,dt,u,varargin{:}); maf=max(mf,maf);
%     [t,x] = ode23(fxu,tv,xi,odeset('MaxStep',dt,'InitialStep',dt),u,0,varargin);   %V5.1
%     keyboard
%     [t,x] = ode23(fxu,tv,xi);       %V4.2c2
     [n,m]=size(x); xi=x(n,:)'; xk=[xk; xi'];
     xout=[xout; x(1:n-1,:)]; tout=[tout; tv(1:n-1)'];
     if flgy;
       yout=[yout; y(1:n-1,:)];
       if i==1; yk=[y(1,:); y(n,:)];
       else yk=[yk; y(n,:)]; end
     end
     if fileflag('digsim.txt');
       if df; close(hw); end;
       error(' Termination due to digsim.txt');
     end;
     if df; waitbar(i/(nt-1)); end;
  end
  if df; close(hw); end
  xout=[xout; xi']; yout=[yout; y(end,:)]; tout=[tout; tv(n)];