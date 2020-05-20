function simopt(sysdyn,sysout,phi,fname,optmat,uydat)
% SIMOPT : Function that simulates the digital optimal control and computes
%          the corresponding output
%
%          simopt(sysdyn,sysout,phi,fname,optmat,uydat)
%
%          Input:
%          sysdyn : String with name of function containing the system dynamics
%                   format function [dxdt]=f(t,x,u,varargin)
%          sysout : String with name of function containing the system dynamics
%                   format function [y]=g(t,x,u,varargin)
%          phi    : String containing the name of phi the terminal costs
%                   format function [phi]=phi(t,x(tf))
%          fname  : String containing the name of the file that contains
%                   the optimal control data
%                   If there are external inputs: extinp=[t d]
%
%          Output
%          optmat : Name of the file that contains the optimal control simulation data
%          uydat  : Name of the ascii file that contains [u y] in the form of a table
%
% GvW 08-02-2010

% Check input
if nargin~=6; error(' 5 input required'); end
if ~(ischar(sysdyn) || isa(sysdyn,'function_handle')); 
  error(' 1st argument must be a string or function handle'); 
end
if ~(ischar(sysout) || isa(sysout,'function_handle')); 
  error(' 2nd argument must be a string or function handle'); 
end
if ~(ischar(phi) || isa(phi,'function_handle')); 
  error(' 3rd argument must be a string or function handle'); 
end
if ~isstr(fname); error(' 4th input must be a string'); end
if ~isstr(optmat); error(' 5th input must be a string'); end
if ~isstr(uydat); error(' 6th input must be a string'); end

% load data & determine dimensions
% !!!!!! varin={[pv] [extinp]} !!!!!!
eval(['load ' fname]); if ~exist('pv','var'); pv=[]; end
if ~exist('extinp','var');
  extinp=[]; dimd=0;
else
  dimd=max(0,size(extinp,2)-1);
end
dt=f_step; tk=tu(:,1); ustar=tu(:,2:end); x0=tx(1,2:end)';
dimxd=size(x0,1)-1; dimu=size(ustar,2);
dimy=length(feval(sysout,0,x0(1:dimxd))); dimx=dimxd-dimd;
xk=tx(:,2:end); lk=tla(:,2:end-1); ub=u_constr;

% Time and piecewise constant control generation
[tu,u]=stair(tk,ustar(1:size(ustar,1)-1,:)); dt=tk(2)-tk(1);

% Plot optimal piecewise constant control
figure; plot(tu,u);
title('Optimal control'); xlabel('time'); ylabel('u(t)');
pause(3);
figure;plot(tk,xk);
title(['Optimal state trajectory, Cost:' num2str(J_optim)]); xlabel('time'); ylabel('x(tk)');
figure;plot(tk,lk);
title(['Optimal costate trajectory, Cost:' num2str(J_optim)]); xlabel('time'); ylabel('l(tk)');
pause(3);
  
% Determine dimensions
dims=[dimx,dimu,dimy,dimd,length(pv)];

% Integration
disp('  Simulation in progress')
[xstar,ystar,tout,xout,yout,maf]=digsim(sysdyn,dims,tk,ustar,x0,dt,sysout,1,pv,extinp);
disp('  Simulation ready');

% Determination final state (xf),
% and minimum costs (cost) 
xf=xstar(end,:)'; cost=xf(dimxd+1)+feval(phi,tout(end),xf(1:dimx))

% Plot of optimal system reponse
figure; plot(tout,xout(:,1:dimx));
title(['Optimal state trajectory, Costs:' num2str(cost)])
xlabel('time [secs]'); ylabel('x(t)');

% Plot of running costs
figure; plot(tout,xout(:,dimx+1));
title(['Optimal running costs, Costs:' num2str(cost)])
xlabel('time [secs]'); ylabel('Jr(t)');
% Determine scale to one factors x,u,y: ds
% Indices of states that are entirely zero : ix0
[dsx,ix0]=sc21mcf(xout); [ds]=[dsx sc21mcf(u) sc21mcf(yout)]; ds=ds';
%ds=[1./max(abs(xout)) 1./max(abs(u)) 1./max(abs(yout))]';

% External inputs handling & scaling factor handling
if dimd;
  % Determine scaling to 1 factors external inputs : dsd
  [dsd]=sc21mcf(extinp(:,2:end)); dsd=dsd'; %
    
  % If their number equals the number of external inputs
  % associate these states to the external inputs. They become
  % deviation variables of the external inputs that may be corrupted
  % with noise. The scaling factors of the external inputs are applied
  % to these states
  if sum(ix0)==dimd;
    ds(ix0)=dsd; ix0=ix0.*[1:length(ix0)]; ix0=ix0(ix0~=0);
    disp(['  States ' num2str(ix0) ' are zero and associated with the external inputs']);
  else
    disp('  No states associated with the external inputs');
  end
  % Compute perturbations for finite differences
  % (not used/needed in case of automatic differentiation
else
  dsd=[]; disp('  No external inputs');
end
% Compute perturbations for finite differences
% (not used/needed in case of automatic differentiation
xudev=1e-3./ds(1:dimxd+dimu+1);

% Save result
eval(['save ' optmat ' dims tk ustar xstar lk ub ystar x0 cost dt tout xout pv ds dsd xudev maf extinp']);
uy=[ustar ystar]; eval(['save -ascii ' uydat ' uy']);