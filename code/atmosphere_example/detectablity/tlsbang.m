% TLSBANG: Optimal bang-bang control with constraints
%          & temporal linear system structure
%
% Input:
% function minimization
%        opt(1) : 1=display intermediate results
%        opt(13): number of equality constraints
%        opt(14): 500, maximum number of function evaluations
%        tguess : array of switching times
%        swi    : associated array of indices of the switching input variables
%        ub     : [ul uu], lower and upperbounds of u
%        ui     : initial values u, elements of ul, ub
%        x0     : initial state system
%        tf     : fimal time
%        p      : parameter vector
%        td     : external inputs [t d], t: columnvector d: stacked row vectors
%
% GvW 28-8-2007


% Problem aand optimization parameters
% 1 System parameters
td=[]; x0=[-10; 10; 0];
tf=4; dt=0.01; N=round(tf/dt)+1; dt=tf/N; % parameters, initial state, final time

% 2 Optimization settings
if p==0; tguess=[1.4109 3.7372]; else tguess=[0.9 3.8]; end % guess of switching times
swi=[1 1]; % associated indices of control input that switches at tguess
ub=[-10 10]; ui=ub(1); % bounds [ul uu] and initial value control
% time lower and upper bounds
tlb=[1e-6*tf]*ones(size(tguess)); tub=[(1-1e-6)*tf]*ones(size(tguess));
opt=[]; opt(1)=1; opt(13)=0; opt(14)=500; % settings for constr

% Dimension checks
if length(tlb)~=length(tguess); error(' tguess and tlb incompatible'); end;
if length(tub)~=length(tguess); error(' tguess and tub incompatible'); end;

% Perform the constraint function minimization and return the solution x
tsw=constr('tlssys',tguess,opt,tlb,tub,[],swi,ub,ui,x0,tf,p,td);

% Obtain the optimal control and optimal state
[f,g,tt,xt,tsw,us,xind]=tlssys(tsw,swi,ub,ui,x0,tf,p,td); f
% Plot the optimal state
hf=figure; plot(tt,xt(:,1:end-1)');
% xlabel('Time'); ylabel('states'); title('Optimal x(t)');
% Plot the optimal control
hold on; stair([0 tsw tf],us);
axis([0 tf min(ub(:,1))-1 max(ub(:,2))+2]);
xlabel('Time'); ylabel('Optimal state & control');
%title('Optimal state & control trajectory');
hold off;
print(hf,'Fig1.emf','-dmeta'); % Save figure
xind=cumsum(xind);

% Compute the costate
% st=0.25; Nc=round(tf/st); st=tf/Nc; ni=20; dt=st/ni; At=zeros(nx,nx*(ni+1)); tim=0;
% ti=[tt(1):(tt(end)-tt(1))/Nc*ni):tt(end)]; xi=interp1(tt,xt,ti,ms); ih=1;
% for i1=1:Nc
%   ip1=1; ip2=nx;
%   for i2=1:ni
%     xc=xi(ih,:)'; x1=xc(1); x2=xc(2); At(:,ip1:ip2)=eval(A);    
%     ih=ih+1; ip1=ip1+nx; ip2=ip2+nx;
%   end
%   [PHIk,GAMk,Qk,Rk,Mk,Vk,gk]=edortv(a,b,q,r,v,ts,z)  
% end