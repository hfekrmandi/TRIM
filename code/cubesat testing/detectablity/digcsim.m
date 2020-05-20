  function [uk,xk,yk,tout,xout]=digcsim(fxu,dflg,tk,uk,yk,x0,dx0,dt,gxu,xfkl,ub,cs,vt,wt,varargin)
%
% DIGCSIM : Function simulates a digital optimal control system
%
%          [uk,xk,yk,tout,xout]=digcsim(fxu,tk,uk,yk,x0,dx0,dt,gxu,xfkl,ub,cs,vt,wt)
%
%          Input:
%                 fxu   :  string with function name xdot = fxu(t,x,flg,u) V5.1
%                          or xdot = fxu(t,x) V4.2c2
%                 dflg  :  display flag if dflg~=0 display
%                 tk    :  column vector with sampling instants
%                 uk    :  piecewise constant control
%                 yk    :  expected outputs
%                 x0    :  initial state
%                 dx0   :  deviations initial state
%                 dt    :  optional prescribed time steps
%                          within sampling interval that determines tout,xout
%                 gxu   :  optional string with output function name
%                          y = gxu(t,x,flg,u)
%                 xfkl  :  Digital compensator
%                 ub    :  Control bounds ub=[umin umax];
%                 cs    :  control correction multiplication factor
%                 vt    :  realization system noise
%                 wt    :  realization measurement noise
%           varargin    :  Additional inputs. Varargin{1:2}={[pv] [extinp]}
%
%	   Output:
%                 uk    :  controls at the sampling instants
%                 xk    :  states at the sampling instants
%                 yk    :  outputs at the sampling instants
%                 tout  :  t0:dt:tN
%                 xout  :  states at tout
%

% check inputs
  if nargin<12; error('  At least 12 input argumenst required'); end;
  if ~(ischar(fxu) || isa(fxu,'function_handle')); 
    error(' 1st argument must be a string or function handle'); 
  end
  [nt,mt]=size(tk); [nu,mu]=size(uk);
  if mt~=1; error('  2nd input must be a column vector'); end;
  if  nu<nt-1; error('  tk and uk incompatible'); end;
  m=size(uk,2);
  if nargin<10; cb=0; ub=[]; end;
  if ~isempty(ub)
    cb=1;
    [mh,nh]=size(ub);
    if mh~=m | nh~=2; error(' ub should be m x 2'); end;
    umin=ub(:,1); umax=ub(:,2);
    if any(umax-umin<=0); error(' Lower bound control exceeds upperbound'); end;
  else
    cb=0;
  end
%  global u  %V4.2c2

% integrate the system and
% save the states at the sampling instants
  xi=x0+dx0; xk=xi'; tout=[]; xout=[]; yout=[]; u=uk(1,:)';
  [nc,ncn,nh,nh,f,k,l,xc]=getfkl(xfkl,1); flg=0; it=0; 
  flag=0; save digcsim.txt flag -ascii;
  if dflg;
    hw=waitbar(0,'Simulation digital control system ...');
  end
  for i=1:nt-1
    % Update times
    ti=tk(i); tf=tk(i+1); th=(tf-ti)/2;
    ni=round((tf-ti)/dt); dt=(tf-ti)/ni;
    % Update output
    [y] = feval(gxu,ti,xi,flg,u,varargin{1:2});
    y=y+wt(:,i); dy=y-yk(i,:)';
    % Update compensator
    du=-l*xc;xc=f*xc+k*dy; u=uk(i,:)'+cs*du;
    if cb; u=min(u,umax); u=max(u,umin); end;
    if i~=nt-1
      [nc,ncn,nh,nh,f,k,l]=getfkl(xfkl,i+1);
    end
    % Update output and control matrices
    yk(i,:)=y'; uk(i,:)=u';
    % Simulate over one sampling interval
    if dt<=0; dt=th; else; min(dt,th); end; tv=ti:dt:tf;
    
    if size(vt,2)==nt-1;
    % Discrete-time system noise addition
      [t,x,y]=ode23fs(fxu,gxu,tv,xi,dt,u,varargin{:});
      [n,m]=size(x); dx=vt(:,i)*[0:1/(size(x,1)-1):1];
      x=x+dx'; xi=x(n,:)'; xk=[xk; xi'];
    else
    % Continuous-time system noise addition      
      j=round((tk(i+1)-tk(i))/dt); vh=vt(:,it+1:it+j); it=it+j;
      [t,x,y]=ode23fsw(fxu,gxu,tv,xi,dt,vh,u,varargin{:});
      [n,m]=size(x); xi=x(n,:)'; xk=[xk; xi'];
    end
    xout=[xout; x(1:n-1,:)]; yout=[yout; y(1:n-1,:)]; tout=[tout; tv(1:n-1)'];
    if fileflag('digcsim.txt');
      if dflg; close(hw); end;
      error(' Termination due to digcsim.txt');
    end;
    if dflg; waitbar(i/(nt-1)); end;
  end
  if dflg; close(hw); end
  xout=[xout; xi']; yout=[yout; y(end,:)]; tout=[tout; tv(n)];