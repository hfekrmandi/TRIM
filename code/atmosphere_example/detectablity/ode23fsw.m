function [t,x,y]=ode23fsw(fxu,gxu,tv,xi,dt,vt,varargin)
% ODE23FSW: Runge Kutta third order numerical integration with fixed step size
%           and white noise disturbances
%
%           [t,x]=ode23fsw(fxu,gxu,tv,xi,dt,vt,varargin)
%         
%           Input:
%           fxu   : Name of the function f(x,u) which describes
%                  the system dynamics (dx/dt=f(x,u)) in state=space form
%           gxu   : Name of function g(x,u)
%           tv    : vector that contains the initial and final time: [ti tf]
%           xi    : Initial state x(ti)
%           dt    : Fixed stepsize numerical integration (rounded internally towards
%                   a divisor of tf-ti)
%           vt    : White noise intensity matrix at each time
%     varargin    : Additional inputs. Varargin{1:3}={[u] [pv] [extinp]}
%
%           Output:
%           t     : time vector [ti ti+dt ti+2*dt ... tf]
%           x     : state matrix [x'(ti); x'(ti+dt); ...; x'(tf)]
%           y     : output matrix [y'(ti); y'(ti+dt); ...; y'(tf)]
%
% Gerard van Willigenburg 3-9-2003
%

if ~isempty(gxu); flgy=1; else flgy=0; end

ti=tv(1); tf=tv(end); ni=round((tf-ti)/dt); dt=(tf-ti)/ni; sqdt=sqrt(dt);
flg=0; lx=length(xi); x=zeros(ni+1,lx); t=zeros(ni+1,1); x(1,:)=xi'; t(1)=ti;
if flgy;
  yi=feval(gxu,ti,xi,varargin{:}); ly=length(yi); y=zeros(ni+1,ly); y(1,:)=yi';
end

for i=1:ni
  xh=xi; th=ti;
  [dx] = feval(fxu,th,xh,varargin{:}); k0=dt*dx;
  if length(dx) ~= lx; error('  Dimensions x and dx do not match'); end;
  xh=xi+0.5*k0; th=ti+0.5*dt;
  [dx] = feval(fxu,th,xh,varargin{:}); k1=dt*dx;
  xh=xi+0.5*k1; th=ti+0.5*dt;
  [dx] = feval(fxu,th,xh,varargin{:}); k2=dt*dx;
  xh=xi+k2; th=ti+dt;
  [dx] = feval(fxu,th,xh,varargin{:}); k3=dt*dx;
  xi=xi+(k0+2*k1+2*k2+k3)/6; ti=ti+dt;
% Add noise to state
  xi=xi+sqdt*vt(:,i);
  x(i+1,:)=xi'; t(i+1)=ti;
  if flgy; yi=feval(gxu,ti,xi,varargin{:}); y(i+1,:)=yi'; end;
end