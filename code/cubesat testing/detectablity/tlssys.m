% TLSSYS: Function to be minimized to obtain the optimal bang-bang control of the system
%
%        Input:
%        ts : array of switching times
%        swi: associated array of indices of the switching input variables
%        ub : [ul ub], lower and upperbounds of u
%        ui : initial values u, elements of ul, ub
%        x0 : initial state system
%        tf : fimal time
%        p  : parameter vector
%        td : external inputs [t d], t: columnvector d: stacked row vectors
%
%        Output:
%        f  : costfunction
%        g  : equality/inequality constraints
%        tt : time column vector
%        xt : associated state response
%        ts : switching times
%        us : control history

function [f,g,tt,xt,ts,us,xind]=tlssys(ts,swi,ub,ui,x0,tf,p,td)

% Dimension checks
if length(ts)~=length(swi);
  error(' ts and swi must have same length');
end
[n1,m1]=size(ui); [n2,m2]=size(ub);
if n1~=m1; error(' ui, ub incompatible'); end;
if m1~=1; error('ui must be columnvector'); end
if m2~=2; error(' ub must be mx2 matrix'); end;

tt=[]; % vector to store time instants obtained from numerical integration
xt=[]; % matrix to store states obtained from numerical integration
xc=[]; % matrix to store states for evaluation of constraints

tk=0; K=2; [ts,ind]=sort(ts); swi=swi(ind);
th=[ts tf]; us=[]; u=ui; xind=[];
% initial time, number of time points, initial state
for k=1:length(th)
  us=[us u]; tk1=th(k); % next time point
  if tk1-tk < 1e-6; tk1=tk+1e-6; end;
  [t,x]=ode45(@tlsdyn,[tk tk1],x0,[],u); % numerical integration
  if k==1; 
    tt=[tt; t]; xt=[xt; x]; % store time points and states
    xind=[xind size(x,1)];
  else
    tt=[tt; t(2:end)]; xt=[xt; x(2:end,:)]; % store time points and states
    xind=[xind size(x,1)-1];
  end
  tk=tk1; x0=x(end,:)'; % update time and the initial state 
  xc=[xc x0]; % update state matrix for constraint evaluation
  if k~=length(th)
    if u(swi(k))==ub(swi(k),1);
      u(swi(k))=ub(swi(k),2);
    else
      u(swi(k))=ub(swi(k),1);
    end;% Switch control
  end
end;
f=x0(3)+(x0(1)-5)*(x0(1)-5)+(x0(2)-5)*(x0(2)-5); g=[];
%g=[x0(1)-5; x0(2)-5]; %Compute J and constraints