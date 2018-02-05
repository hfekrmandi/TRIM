% SIMCSYS : Simulate the digital optimal control system
%
%          L.G. van Willigenburg 31-1-'96.
%
% Input: rp,cw,sd,mulnoi (repeat, noise discrete-time)
%   rp~=0: repeat simulation and switch cs 1 <-> 0 (default 0)
%   cs: control correction multiplication factor (default 1)
%   cw~=0: continuous white system noise generation else discrete-time (default 0)
%   sd : seed random number generator (default 1)
%   mulnoi : Noise multiplication factor (default 1)
%   ix : indices of the state x that are plotted (default 1:dimxd)
%   dflg : display flag, if dflg~=0 display (default 0)

% Load optimal control and compendator data for simulation  
  load uystar; load digcomp;
  dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;

% Default values input parameters
  if ~exist('dflg','var'); dflg=0; end;
  if ~exist('rp','var'); rp=0; end;
  if ~exist('cw','var'); cw=1; end;
  if ~exist('sd','var'); sd=1; end;
  if ~exist('mulnoi','var'); mulnoi=1; end;
  if ~exist('ix','var'); ix=1:dimxd; end

% Initialization seed
  randn('state',sd);
  
% If rp==1 repeat the simulation and toggle the compensation
  if rp && exist('dx0','var');
    if ~exist('cs','var'); cs=1; elseif cs==0; cs=1; else cs=0; end;
    if dflg
      if cs==0; disp(' Repeated simulation without compensation');
      else disp(' Repeated simulation with compensation'); end
    end
  else
%   Generate system, measurement and initial state noise
    cs=1; if dflg; disp(' New simulation with compensation'); end;
    [x0,dx0,vh,wh]=gennoi(dims,lqgarr,tk,dt,x0,vt,wt,cw,mulnoi);
  end
 
% Integration
%  disp('  Simulation in progress')
  [umk,xmk,ymk,tmout,xmout]=digcsim(sysdyn,dflg,tk,ustar,ystar,x0,dx0,dt,sysout,xfkl,1.2*ub,cs,vh,wh,pv,extinp,1);
%  disp('  Simulation ready');
  
% Determination final state (xf),
% and minimum costs (cost) 
  xf=xmk(end,1:dimxd+1)'; costr=xf(end)+feval(phi,tk(end),xf(1:dimxd));

% Time and piecewise constant control generation
  [tu,um]=stair(tk,umk(1:end-1,:));
  [tu,us]=stair(tk,ustar(1:end-1,:)); 

  if dflg;
% Plot optimal and real piecewise constant control
    figure(1); plot(tu,um*diag(ds(dimxd+2:dimxd+dimu+1)),tu,us*diag(ds(dimxd+2:dimxd+dimu+1)));
    title('Real and optimal scaled control'); xlabel('time [secs]'); ylabel('u(t)');
    pause(3);

% Plot of optimal and real system response
    figure(gcf+1); plot(tmout,xmout(:,ix)*diag(ds(ix)),tout,xout(:,ix)*diag(ds(ix)));
    title(['Real scaled state trajectory, Real Costs:' num2str(costr)])
    xlabel('time [secs]'); ylabel('x(t)');

% Plot of optimal and real running costs  
    figure(gcf+1); plot(tmout,xmout(:,dimxd+1),tout,xout(:,dimxd+1));
    title(['Real and optimal running costs, Costs:' num2str(costr)])
    xlabel('time [secs]'); ylabel('Jr(t)');
    disp([' Real cost : ' num2str(costr)]);
  end