function [tu,tx,tla,tt,xt,yt,J_optim,x0,u_constr,probname,varargin]=opconstr1(probname,startname,optname,fopt,m,varargin)
% OPCONSTR1: Solve a digital optimal control problem with fixed final time
%            and bounded control using CONSTR and the user supplied gradient dH/du.
%
%            [u,s,tsar]=opconstr(probname,startname,optname,fopt,m,varargin) or
%            [u,s,tsar]=opconstr(probname,tu       ,optname,fopt,m,varargin)
%
% Input:
% probname  : name of file containing the
%             digital optimal control problem data
% startname : Name of mat-file containing an initial guess.
%             If the file does not exist or [] is specified
%             default initial values are used.
% optname   : Name of mat-file with solution i.e. the a,b,l and T values.
% fopt      : Options for function constr.m
% m         : The figure number for display of the solution (if empty no display)
% varargin  : Additional parameters passed to the problem file
%             varargin{1:3}={[tu] [pv] [extinp]}
%
% Output:
% tu        : time and optimal control matrix [[t0;t1..;tN],[u0';u1';..uN']];
% tx        : time and optimal state matrix [[t0;t1..;tN],[x0';x1';..xN']];
% tla       : time and optimal costate matrix [[t0;t1..;tN],[l0';l1';..lN']];
% tt        : Continuous time
% xt        : Intersample states
% yt        : Intersample outputs
% J_optim   : optimal cost
% x0        : Initial state
% u_constr  : Control bounds
% probname  : problem file name
%
%            Edit breakdop.txt into non-zero to terminate execution
%            Ctrl-C, Ctrl-Break give problems
%
%            GvW Last modified : 4-6-2001

global J_glob tol_glob f_step; J_glob=1e60; tol_glob=1e60; f_step=0;

% Check inputs
if nargin<3; error('At least 3 input arguments required');
elseif nargin==3; fopt=zeros(5,1); m=[];
elseif nargin==4; m=[]; end
if isempty(optname)&nargout==0
  error('No output file and no output arguments specified');
end
if ~isempty(startname);
  if ~isstr(startname);
    tu=startname; startname='';
  end;
end
if ~isempty(optname);
  if ~isstr(optname);
    error('3rd argument must be a string');
  end;
end
if length(m)==0;
  m=0;
elseif length(m)~=1;
  error('5th input argument must be an integer > 1');
else
  if m<1 || m-round(m)~=0
    error('5th input argument must be an integer > 1');
  end
end
fopt=[fopt(:); zeros(5-length(fopt),1)];

if exist([startname '.mat'])
  load(startname);
  if isempty(varargin) && exist('varin','var')
      varargin=varin;
  end;
elseif ~exist('tu','var') || isempty(tu)
  % Get initial control (flag=4)
  [tu,f_step]=feval(probname,[],[],[],[],4,varargin{:});
end
if f_step==0
  [h1,f_step]=feval(probname,[],[],[],[],4,varargin{:});
end

% Get initial state (flag=5)
[x0]=feval(probname,[],[],[],[],5,varargin{:});

% Get system parameter vector (flag=6)
[u_constr,du_max]=feval(probname,[],[],[],[],6,varargin{:});

if isempty(u_constr);
  if isempty(du_max)
     error(' Either control bounds or a maximum control change must be specified')
  else
    u_constr=[-10*du_max; 10*du_max]; disp(' OPCONSTR : Control bounds set to +/-10 times maximum control change');
  end
end

% Create control bound for constr.m
ulow=[]; uup=[];
for i=1:size(tu,1)-1
  ulow=[ulow u_constr(:,1)']; uup=[uup u_constr(:,2)'];
end

% Get number of equality constraints in PSI (flag=8)
[neq]=feval(probname,[],[],[],[],8,varargin{:});
if isempty(neq); error(' Number of equality constraints in PSI not specified'); end;
fopt(13)=neq; ctol=fopt(4);

% Check if computation should be performed without gradients
[nograd]=feval(probname,[],[],[],[],9,varargin{:});
if isempty(nograd); nograd=0; end;
if nograd; grstr=[]; else grstr='grdHdu'; end

% Optimization
tic;
[u,fopt]=constr('minju',tu(1:end-1,2:end)',fopt,ulow,uup,grstr,probname,x0,tu(:,1)',u_constr,neq,ctol,varargin{:});
toc; if fopt(19)~=1; load utemp; end

% Compute costs, states, and dH/du of solution
[J_optim,G,x]=minju(u,probname,x0,tu(:,1)',u_constr,neq,ctol,varargin{:});
[Hu,HuPSI,lat]=grdHdu(u,probname,x0,tu(:,1)',u_constr,neq,ctol,varargin{:});

% Adapt to output format 
tu=[tu(:,1),[u'; zeros(1,size(tu,2)-1)]]; tx=[tu(:,1),x']; tla=[tu(:,1),lat'];
  
if ~isempty(optname)
  varin=varargin;
  save(optname,'J_optim','x0','tu','tx','tla','u_constr','f_step','varin','probname');
  disp(['Solution saved in ' optname '.mat']);
end

% Simulation optimal solution to obtain intersample behavior
[tt,xt,yt]=digsimo(u,probname,x0,tu(:,1)',varargin{:});

if m~=0
  % Plot state trajectory and output
  figure(m); hold off; clf; plot(tt,xt(:,1:2));
  hold on;  stair(tu(:,1)',tu(1:end-1,2:end)','m'); hold off;
  xlabel('Time'); ylabel('States');
  str=['Optimal control & state trajectory min(J)=' num2str(J_optim)];
  str=[str ' tf=' num2str(tu(end,1))]; title(str);

  figure(m+1); hold off; clf; plot(tt,yt);
  xlabel('Time'); ylabel('Outputs');
  str=['Optimal output trajectory min(J)=' num2str(J_optim)];
  str=[str ' tf=' num2str(tu(end,1))]; title(str);

  % Plot digital optimal control
  figure(m+2); hold off; clf; stair(tu(:,1)',tu(1:end-1,2:end)');
  xlabel('Time'); ylabel('Control');
  str=['Optimal control min(J)=' num2str(J_optim)];
  str=[str ' tf=' num2str(tu(end,1))]; title(str);
end