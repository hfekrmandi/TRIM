function [A,B]=dlfxuad(fcm,t,ts,x,dt,u,varargin)
% DLFXUAD : Function differentiates the integrated vector function fcn(x,u)
%           over t to t+ts with respect to x, result A, and u, result B.
%           A,B is the linearised equivalent discrete-time system over [t,t+ts].
%           Automatic differentiation is employed.
%
%           [A,B]=dlfxuad(fcm,t,ts,x,dt,u,varargin)
%
%           Input   :
%           fcm     : f(x,u), m-function string
%           t       : begin sampling period
%           ts      : length sampling interval
%           x       : state
%           dt      : numerical integration time step
%           u       : control input
%    varargin       : additional inputs. Varargin{1:2}={[pv] [extinp]}
%
%           Output  :
%           A,B     : Linearised equivalent discrete-time system over [t,t+ts]
%
% GvW : Last updated 4-6-2001

%
if ~isstr(fcm); error('  fcm must be a string with a function name'); end;

[nx,nh]=size(t);
if nx~=1 || nh~=1; error('  t must be a scalar'); end;

[nx,nh]=size(ts);
if nx~=1 || nh~=1; error('  ts must be a scalar'); end;

[nx,nh]=size(x);
if nh~=1; error('  x must be a column vector'); end;

[nt,nh]=size(t);
if nt~=1 || nh~=1; error('  dt must be a scalar'); end;

[nu,nh]=size(u);
if nh~=1; error('  u must be a column vector'); end;

% Gradient computation (intlab)
if exist('gradientinit','file') 
  % Compute trajectory x(t) over ts
  [tim,xt] = ode23fs(fcm,[],[t t+ts],x,dt,u,varargin{:});

  % Get A=df/dx,  B=df/du over ts
  A=eye(nx); B=zeros(nx,nu);
  e=A; nt=size(xt,1); 
  delt=ts/(nt-1); delta=0.5*delt;
  ix=1:nx; iu=nx+1:nx+nu;

  for i=1:nt
    % Evaluate A(t)
    th=t+(i-1)/(nt-1)*ts; xu=gradientinit([xt(i,:)';u]);
    dxts=feval(fcm,th,xu(ix),xu(iu),varargin{:});
    %dxts=dynex(th,xh,u,pv);
    AB=dxts.dx; Atn=AB(:,ix); Bt=AB(:,iu);
    % Update A, B
    if i>1
      ad=delta*(Atn+At); Aud=e+ad+0.5*ad*ad; A=Aud*A;
      B=Aud*(B+Bud); Bud=delta*Bt; B=B+Bud;
    else
      Bud=delta*Bt;
    end
    % Replace old value At by new one
    At=Atn;
  end
   
% Gradient computation (tomlab)
elseif exist('fmad','file');
  % Compute trajectory x(t) over ts
  [tim,xt] = ode23fs(fcm,[],[t t+ts],x,dt,u,varargin{:});  

  % Get A=df/dx,  B=df/du over ts
  A=eye(nx); B=zeros(nx,nu);
  e=A; nt=size(xt,1); 
  delt=ts/(nt-1); delta=0.5*delt;
  ix=1:nx; iu=nx+1:nx+nu;
  
  for i=1:nt
    % Evaluate A(t)
    th=t+(i-1)/(nt-1)*ts; xu=fmad([xt(i,:)';u],eye(nx+nu));
    dxts=feval(fcm,th,xu(ix),xu(iu),varargin{:});
    %dxts=dynex(th,xh,u,pv);
    AB=full(getinternalderivs(dxts));
    Atn=AB(:,ix); Bt=AB(:,iu);
    % Update A, B
    if i>1
      ad=delta*(Atn+At); Aud=e+ad+0.5*ad*ad; A=Aud*A;
      B=Aud*(B+Bud); Bud=delta*Bt; B=B+Bud;
    else
      Bud=delta*Bt;
    end
    % Replace old value At by new one
    At=Atn;
  end
  
% No automatic differentiation available use finite differences
else
  % Compute trajectory x(t) over ts
  [tim,xt] = ode23fs(fcm,[],[t t+ts],x,dt,u,varargin{:});  

  % Get A=df/dx,  B=df/du over ts
  A=eye(nx); B=zeros(nx,nu);
  e=A; nt=size(xt,1); 
  delt=ts/(nt-1); delta=0.5*delt;
  ix=1:nx; iu=nx+1:nx+nu;
  
  for i=1:nt
    % Evaluate A(t)
    th=t+(i-1)/(nt-1)*ts; 
    %[Atn,Bt]=linvfxu(fcm,th,xt(i,:)',u,[],[],[],[],[],varargin{:});
    [Atn,Bt]=linfxu1(fcm,th,xt(i,:)',u,[],[],[],varargin{:});
    % Update A, B
    if i>1
      ad=delta*(Atn+At); Aud=e+ad+0.5*ad*ad; A=Aud*A;
      B=Aud*(B+Bud); Bud=delta*Bt; B=B+Bud;
    else
      Bud=delta*Bt;
    end
    % Replace old value At by new one
    At=Atn;
  end
end
