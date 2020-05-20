function [x,OPTIONS,lambda,HESS]=nlconstold(FUNfcn,x,OPTIONS,VLB,VUB,GRADfcn,varargin)
%NLCONST Helper function to find the constrained minimum of a function 
%   of several variables. Called by CONSTR, ATTGOAL. SEMINF and MINIMAX.
%
%   [X,OPTIONS,LAMBDA,HESS]=NLCONST('FUN',X0,OPTIONS,VLB,VUB,'GRADFUN',...
%   varargin{:}) starts at X0 and finds a constrained minimum to 
%   the function which is described in FUN. FUN is a four element cell array
%   set up by PREFCNCHK.  It contains the call to the objective/constraint
%   function, the gradients of the objective/constraint functions, the
%   calling type (used by OPTEVAL), and the calling function name. 
%
%   Copyright 1990-2002 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2002/03/12 20:36:14 $
%   Andy Grace 7-9-90, Mary Ann Branch 9-30-96.
%
%   $Revised by GvW 2-5-2007:
%   OPTIONS(19)=1 indicates succesful termination.
%   See also modified constopt.m.

%   Called by CONSTR, SEMINF, ATTGOAL, MINIMAX.
%   Calls OPTEVAL.
%

% Expectations: GRADfcn must be [] if it does not exist.
global OPT_STOP OPT_STEP;

% GvW part
OPTIONS(19)=0; % OPTIONS(19)=1 if termination succesful default 0 GvW
nfb=0; IHNU=0; Hesstop=0; fs1=0; fs2=0; % counters
if length(OPTIONS)>=20 && OPTIONS(20)~=0;
    merfac=OPTIONS(20);
else
    merfac=1;
end    % counter for convergence, merit factor GvW

if length(OPTIONS)>=23 && OPTIONS(23)~=0;
    Hessmax=OPTIONS(23);
else
    Hessmax=5;
end    % counter for Hessian modified twice
% end GvW part

if length(OPTIONS)>=24 && OPTIONS(24)~=0;
    fmax1=OPTIONS(24);
else
    fmax1=10;
end    % counter for improvement f
% end GvW part

if length(OPTIONS)>=25 && OPTIONS(25)~=0;
    fmax2=OPTIONS(25);
else
    fmax2=10;
end    % counter for improvement f
% end GvW part

OPT_STEP = 1; 
OPT_STOP = 0; 
% Initialize so if OPT_STOP these have values
lambda = []; HESS = [];

Nlconst = 'nlconst';

% Set up parameters.
XOUT=x(:);

VLB=VLB(:); lenvlb=length(VLB);
VUB=VUB(:); lenvub=length(VUB);
%bestf = Inf; 
fmem=1e110; bestf=1e100; bestg=bestf;%GvW

nvars = length(XOUT);

CHG = 1e-7*abs(XOUT)+1e-7*ones(nvars,1);
if lenvlb*lenvlb>0
      if any(VLB( (1:lenvub)' ) > VUB), error('Bounds Infeasible'), end
end
for i=1:lenvlb
       if lenvlb>0,if XOUT(i)<VLB(i),XOUT(i)=VLB(i)+1e-4; end,end
end
for i=1:lenvub
       if lenvub>0,if XOUT(i)>VUB(i),XOUT(i)=VUB(i);CHG(i)=-CHG(i);end,end
end

t0=cputime; %GvW
% Used for semi-infinite optimization:
s = nan; POINT =[]; NEWLAMBDA =[]; LAMBDA = []; NPOINT =[]; FLAG = 2;
OLDLAMBDA = [];

x(:) = XOUT;  % Set x to have user expected size
% Compute the objective function and constraints
if strcmp(FUNfcn{4},'seminf')
  [f,g,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
     semifun(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
else
  [f,g,msg] = opteval(x,FUNfcn,varargin{:});
  error(msg);
  g = g(:);
end
f=merfac*f; %GvW
if isempty(f)
  error('FUN must return a non-empty objective function.')
end
ncstr = length(g);

% Evaluate gradients and check size
if isempty(GRADfcn)  
  analytic_gradient = 0;
else
  analytic_gradient = 1;
  [gf_user,gg_user,msg] = opteval(x,GRADfcn,varargin{:});
  error(msg);            
  gf_user = gf_user(:);
  % Both might evaluate to empty when expression syntax is used
  if isempty(gf_user) && isempty(gg_user) 
    analytic_gradient = 0;
  else  % Either gf or gg is defined
     if length(gf_user) ~= nvars
       error('The objective gradient is the wrong size.')
     end
     if isempty(gg_user) && isempty(g)
       % Make gg compatible
       gg = g'; gg_user = g';
     else % Check size of gg
       [ggrow, ggcol] = size(gg_user);
       if ggrow ~= nvars 
          error('The constraint gradient has the wrong number of rows.')
       end
       if ggcol ~= ncstr
          error('The constraint gradient has the wrong number of columns.')
       end
     end % isempty(gg_user)
  end % isempty(gf_user) && isempty(gg_user)
end % isempty(GRADfcn) 
 
OLDX=XOUT;
OLDG=g;
OLDgf=zeros(nvars,1);
gf=zeros(nvars,1);
OLDAN=zeros(ncstr,nvars);
LAMBDA=zeros(ncstr,1);
sizep = length(OPTIONS);
OPTIONS = constopt(OPTIONS);
if lenvlb*lenvlb>0
      if any(VLB((1:lenvub)') > VUB), error('Bounds Infeasible'), end
end
for i=1:lenvlb
       if lenvlb>0,if XOUT(i)<VLB(i),XOUT(i)=VLB(i)+eps; end,end
end
OPTIONS(18)=1;
if OPTIONS(1)>0
   if OPTIONS(7)==1
        disp('')
        disp('f-COUNT     MAX{g}         STEP  Procedures');
   else
    disp('')
        disp('f-COUNT   FUNCTION       MAX{g}         STEP  Procedures');
   end
end
HESS=eye(nvars,nvars);
if sizep<1 ||OPTIONS(14)==0, OPTIONS(14)=nvars*100;end
OPTIONS(10)=1;
OPTIONS(11)=1;
GNEW=1e8*CHG;


%---------------------------------Main Loop-----------------------------
status = 0; 
while status ~= 1

%----------------GRADIENTS----------------

  if ~analytic_gradient || OPTIONS(9)
% Finite Difference gradients (even if just checking analytical)
     POINT = NPOINT; 
     oldf = f;
     oldg = g;
     ncstr = length(g);
     FLAG = 0; % For semi-infinite
     gg = zeros(nvars, ncstr);  % For semi-infinite
% Try to make the finite differences equal to 1e-8.
     CHG = -1e-8./(GNEW+eps);
     CHG = sign(CHG+eps).*min(max(abs(CHG),OPTIONS(16)),OPTIONS(17));
     OPT_STEP = 1;
     for gcnt=1:nvars
        if gcnt == nvars, 
           FLAG = -1; 
        end
        temp = XOUT(gcnt);
        XOUT(gcnt)= temp + CHG(gcnt);
        x(:) =XOUT; 
            
        if strcmp(FUNfcn{4},'seminf')
          [f,g,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
             semifun(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
        else
          [f,g,msg] = opteval(x,FUNfcn,varargin{:});
          error(msg);
          g = g(:);
        end
        f=merfac*f; %GvW
        OPT_STEP = 0;
        if OPT_STOP
           break;
        end
        % Next line used for problems with varying number of constraints
        if ncstr~=length(g), 
           diff=length(g); 
           g=v2sort(oldg,g); 
        end

        gf(gcnt,1) = (f-oldf)/CHG(gcnt);
        if ~isempty(g)
           gg(gcnt,:) = (g - oldg)'/CHG(gcnt); 
        end
        XOUT(gcnt) = temp;
     end % for 
     if OPT_STOP
        break;
     end
          
% Gradient check
     if OPTIONS(9) == 1 && analytic_gradient
        gfFD = gf;
        ggFD = gg; 
        gg = gg_user;
        gf = gf_user;

        disp('Function derivative')
        if isa(GRADfcn{1},'inline')
          graderr(gfFD, gf, formula(GRADfcn{1}));
        else
          graderr(gfFD, gf, GRADfcn{1});
        end
        if ~isempty(gg)
          disp('Constraint derivative')
          if isa(GRADfcn{3},'inline')
            graderr(ggFD, gg, formula(GRADfcn{3}));
          else
            graderr(ggFD, gg, GRADfcn{3});
          end
        end
        OPTIONS(9) = 0;
     end % OPTIONS(9) == 1 && analytic_gradient
     FLAG = 1; % For semi-infinite
     OPTIONS(10) = OPTIONS(10) + nvars;
     f=oldf;
     g=oldg;
  else % analytic_gradient && options(9)=0
     % User-supplied gradients
     % gf and gg already computed first time through loop
     if OPTIONS(10) > 1
       gg = zeros(nvars, ncstr);
       [gf,gg,msg] = opteval(x,GRADfcn,varargin{:});
       error(msg);
       gf = gf(:);
       if isempty(gg) && isempty(g)
         gg = g';
       end
     else
       % First time through loop
       gg = gg_user;
       gf = gf_user;
     end
               
     if OPT_STOP
        break;
     end
          
  end  % if ~analytic_gradient || OPTIONS(9)
  AN=gg';
  how='';
  OPT_STEP = 2;

%-------------SEARCH DIRECTION---------------

  for i=1:OPTIONS(13) 
     schg=AN(i,:)*gf;
     if schg>0
       AN(i,:)=-AN(i,:);
       g(i)=-g(i);
     end
  end

  if OPTIONS(11)>1  % Check for first call    
     % For equality constraints make gradient face in 
     % opposite direction to function gradient.
     if OPTIONS(7)~=5,   
       NEWLAMBDA=LAMBDA; 
     end
     [ma,na] = size(AN);
     GNEW=gf+AN'*NEWLAMBDA;
     GOLD=OLDgf+OLDAN'*LAMBDA;
     YL=GNEW-GOLD;
     sdiff=XOUT-OLDX;
     % Make sure Hessian is positive definite in update.
     if YL'*sdiff<OPTIONS(18)^2*1e-3
        while YL'*sdiff<-1e-5
           [YMAX,YIND]=min(YL.*sdiff);
           YL(YIND)=YL(YIND)/2;
        end
        if YL'*sdiff < (eps*norm(HESS,'fro'));
           how=' Hessian modified twice';
           FACTOR=AN'*g - OLDAN'*OLDG;
           FACTOR=FACTOR.*(sdiff.*FACTOR>0).*(YL.*sdiff<=eps);
           WT=1e-2;
           if max(abs(FACTOR))==0; FACTOR=1e-5*sign(sdiff); end
              while YL'*sdiff < (eps*norm(HESS,'fro')) && WT < 1/eps
                 YL=YL+WT*FACTOR;
                 WT=WT*2;
              end
           else
              how=' Hessian modified';
        end
     end

%----------Perform BFGS Update If YL'S Is Positive---------
     if YL'*sdiff>eps
       HESS=HESS ...
      +(YL*YL')/(YL'*sdiff)-((HESS*sdiff)*(sdiff'*HESS'))/(sdiff'*HESS*sdiff);
       % BFGS Update using Cholesky factorization  of Gill, Murray and Wright.
       % In practice this was less robust than above method and slower. 
       %   R=chol(HESS); 
       %   s2=R*S; y=R'\YL; 
       %   W=eye(nvars,nvars)-(s2'*s2)\(s2*s2') + (y'*s2)\(y*y');
       %   HESS=R'*W*R;
     else
       how=' Hessian not updated';
     end

  else % First call
     OLDLAMBDA=(eps+gf'*gf)*ones(ncstr,1)./(sum(AN'.*AN')'+eps) ;
  end % if OPTIONS(11)>1
  OPTIONS(11)=OPTIONS(11)+1;

  LOLD=LAMBDA;
  OLDAN=AN;
  OLDgf=gf;
  OLDG=g;
  OLDF=f;
  OLDX=XOUT;
  XN=zeros(nvars,1);
  if (OPTIONS(7)>0&&OPTIONS(7)<5)
    % Minimax and attgoal problems have special Hessian:
    HESS(nvars,1:nvars)=zeros(1,nvars);
    HESS(1:nvars,nvars)=zeros(nvars,1);
    HESS(nvars,nvars)=1e-8*norm(HESS,'inf');
    XN(nvars)=max(g); % Make a feasible solution for qp
  end
  if lenvlb>0,
    AN=[AN;-eye(lenvlb,nvars)];
    GT=[g;-XOUT((1:lenvlb)')+VLB];
  else
    GT=g;
  end
  if lenvub>0
     AN=[AN;eye(lenvub,nvars)];
     GT=[GT;XOUT((1:lenvub)')-VUB];
  end
  HESS = (HESS + HESS')*0.5;
  
  try
    warning('off');
    [SD,lambda,howqp] = qpsubold(HESS,gf,AN,-GT,[],[],XN,OPTIONS(13),-1, ...
                                 Nlconst,size(AN,1),nvars);
    warning('on');                         
  catch
    SD=[]; warning('on');
  end
   
%  Added by GvW to overcome failures of qpsubold                         
  iHESS=1;
  while isempty(SD) && iHESS<=5;
    disp([' nlconstold: qpsubold failure ' num2str(iHESS)]);
    HESS=HESS+1e-6*norm(HESS)*psd(size(HESS,1));
    try
      warning('off');
      [SD,lambda,howqp] = qpsubold(HESS,gf,AN,-GT,[],[],XN,OPTIONS(13),-1, ...
                                   Nlconst,size(AN,1),nvars);
      warning('on');                         
    catch
      SD=[]; warning('on');
    end
    iHESS=iHESS+1;
    if isempty(SD) && iHESS==6;
      error(' nlconstold terminated due to 5 consequtive failures qpsubold');
    else
      disp(' nlconstold: continuation')
    end
  end
  
  if strcmp(how,' Hessian not updated'); IHNU=IHNU+1; else; IHNU=0; end
  if IHNU>3; disp(' Optimization terminated: No Hessian update'); status=1; end
  
% Continuation original code   
  lambda((1:OPTIONS(13))') = abs(lambda( (1:OPTIONS(13))' ));
  ga=[abs(g( (1:OPTIONS(13))' )) ; g( (OPTIONS(13)+1:ncstr)' ) ];
  if ~isempty(g)
    mg=max(ga);
  else
    mg = 0;
  end
    
  if OPTIONS(1)>0
     if strncmp(howqp,'ok',2); 
       howqp =''; 
     end
        if ~isempty(how) && ~isempty(howqp) 
           how = [how,'; '];
        end
        if OPTIONS(7)==1,
          gamma = mg+f;
          disp([sprintf('%5.0f %12.6g ',OPTIONS(10),gamma), ...
                sprintf('%12.3g  ',OPTIONS(18)),how, ' ',howqp]);
        else
          disp([sprintf('%5.0f %12.6g %12.6g ',OPTIONS(10),f/merfac,mg), ... %GvW
                sprintf('%12.3g  ',OPTIONS(18)),how, ' ',howqp]);
        end
  end
  LAMBDA=lambda((1:ncstr)');
  OLDLAMBDA=max([LAMBDA';0.5*(LAMBDA+OLDLAMBDA)'])' ;

%---------------LINESEARCH--------------------
  MATX=XOUT;
  MATL = f+sum(OLDLAMBDA.*(ga>0).*ga) + 1e-30;
  infeas = strncmp(howqp,'i',1);
  if OPTIONS(7)==0 || OPTIONS(7) == 5
     % This merit function looks for improvement in either the constraint
     % or the objective function unless the sub-problem is infeasible in which
     % case only a reduction in the maximum constraint is tolerated.
     % This less "stringent" merit function has produced faster convergence in
     % a large number of problems.
     if mg > 0
       MATL2 = mg;
     elseif f >=0 
       MATL2 = -1/(f+1);
     else 
       MATL2 = 0;
     end
     if ~infeas && f < 0
       MATL2 = MATL2 + f - 1;
     end
   else
     % Merit function used for MINIMAX or ATTGOAL problems.
     MATL2=mg+f;
  end
%   if mg < eps && f < bestf
   if (mg < OPTIONS(4) && f < bestf) || (bestg > OPTIONS(4) && mg < bestg) %GvW
      bestf = f; bestg=mg;
      bestx = XOUT;
      bestHess = HESS;
      bestlambda = lambda;
   end
   MERIT = MATL + 1;
   MERIT2 = MATL2 + 1; 
   OPTIONS(18)=2;
   while  (MERIT2 > MATL2) && (MERIT > MATL) ...
          && OPTIONS(10) < OPTIONS(14) && ~OPT_STOP
      OPTIONS(18)=OPTIONS(18)/2;
      if OPTIONS(18) < 1e-4,  
        OPTIONS(18) = -OPTIONS(18); 

        % Semi-infinite may have changing sampling interval
        % so avoid too stringent check for improvement
        if OPTIONS(7) == 5, 
          OPTIONS(18) = -OPTIONS(18); 
          MATL2 = MATL2 + 10; 
        end
      end
      XOUT = MATX + OPTIONS(18)*SD;
      x(:)=XOUT; 

      if strcmp(FUNfcn{4},'seminf')
         [f,g,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
              semifun(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
      else
        [f,g,msg] = opteval(x,FUNfcn,varargin{:});
        error(msg);
        g = g(:);
      end
      f=merfac*f; %GvW
      if OPT_STOP
        break;
      end
  
      OPTIONS(10) = OPTIONS(10) + 1;
      ga=[abs(g( (1:OPTIONS(13))' )) ; g( (OPTIONS(13)+1:length(g))' )];
      if ~isempty(g)
         mg=max(ga);
      else
         mg = 0;
      end

      MERIT = f+sum(OLDLAMBDA.*(ga>0).*ga);
      if OPTIONS(7)==0 || OPTIONS(7) == 5
        if mg > 0
          MERIT2 = mg;
        elseif f >=0 
          MERIT2 = -1/(f+1);
        else 
          MERIT2 = 0;
        end
        if ~infeas && f < 0
          MERIT2 = MERIT2 + f - 1;
        end
      else
        MERIT2=mg+f;
      end
   end
%------------Finished Line Search-------------

   if OPTIONS(7)~=5
      mf=abs(OPTIONS(18));
      LAMBDA=mf*LAMBDA+(1-mf)*LOLD;
   end
   % Test stopping conditions (convergence)
   if max(abs(SD))<2*OPTIONS(2) && abs(gf'*SD)<2*OPTIONS(3)*merfac && ...
          (mg<OPTIONS(4) || (strncmp(howqp,'i',1) && mg > 0 ) )
      if OPTIONS(1)>0
         if OPTIONS(7)==1,
            gamma = mg+f;
            disp([sprintf('%5.0f %12.6g ',OPTIONS(10),gamma),...
                  sprintf('%12.3g  ',OPTIONS(18)),how, ' ',howqp]);
         else
            disp([sprintf('%5.0f %12.6g %12.6g ',OPTIONS(10),f/merfac,mg),... %GvW
                  sprintf('%12.3g  ',OPTIONS(18)),how, ' ',howqp]);
         end
         if ~strncmp(howqp, 'i', 1) 
            disp('Optimization Converged Successfully')
            OPTIONS(19)=1; % GvW 7-5-2007 OPTIONS(19)=1 if convergence successful
% The next code is removed because it is incorporated in maximum g GvW                        
%            active_const = find(LAMBDA>0);
%            if active_const 
%                 disp('Active Constraints:'), % Original code
%                 disp(active_const) % Original code
%            else % active_const == 0
%               disp(' No Active Constraints');
%            end 
         end
      end
      if (strncmp(howqp, 'i',1) && mg > 0)
          disp('Warning: No feasible solution found.')
      end
      status=1;
   elseif OPTIONS(10)-nfb > 0.5*OPTIONS(14) || (OPTIONS(10)-nfb > 0.1*OPTIONS(14) && bestg <= OPTIONS(4))%GvW
      if fmem-bestf~=0 && (fmem-bestf)<OPTIONS(3)*merfac && mg>1e3*OPTIONS(4); %GvW
        status=1; % GvW
        disp('Optimization terminated due to lack of improvement'); %GvW
      end %GvW
      nfb=OPTIONS(10);
      if bestf~=1e100; fmem=bestf; end %GvW
   else  % continue
      % NEED=[LAMBDA>0] || G>0
      if OPTIONS(10) >= OPTIONS(14)  || OPT_STOP
         XOUT = MATX;
         f = OLDF;
         if ~OPT_STOP
            disp('Maximum number of function evaluations exceeded;')
            disp('increase OPTIONS(14)')
         end
         status=1;
      end
   end
   if abs(OPTIONS(18))<1e-4; Hesstop=Hesstop+1; else Hesstop=0; end % GvW
   if Hesstop>Hessmax; Hesstop=-1; status=1; end % GvW
   if f>=bestf && mg<OPTIONS(4); fs1=fs1+1; else fs1=0; end % GvW
   if fs1>fmax1; fs1=-1; status=1; end % GvW
   if f>=bestf && mg>OPTIONS(4) && mg>=bestg; fs2=fs2+1; else fs2=0; end % GvW
   if fs2>fmax2; fs1=-1; status=1; end % GvW
end % while status ~= 1

% Update 
OPTIONS(12) = OPTIONS(11);

% If a better unconstrained solution was found earlier, use it:
% if f > bestf % Original
if f > bestf || OPTIONS(19)~= 1 || mg > OPTIONS(4) %GvW
    if bestf~=1e100 %GvW
      XOUT = bestx;
      f = bestf;
      mg = bestg; g=mg; %GvW
      HESS = bestHess;
      lambda = bestlambda;
    end %GvW
end
f=f/merfac; %GvW
OPTIONS(8)=f;
x(:) = XOUT;
if (OPT_STOP)
    disp('Optimization terminated prematurely by user');
end
if mg > OPTIONS(4) % GvW
    disp('No solution found that satisfies the constraints'); % GvW
end % GvW
if Hesstop==-1; % GvW
    disp('Optimization terminated due to repeated small stepsize'); % GvW
end %GvW
if fs1==-1 % GvW
    disp('Optimization terminated due to repeated lack of improvement'); % GvW
end %GvW
tf=cputime; OPTIONS(21)=tf-t0; OPTIONS(22)=mg; %GvW
disp(['CPU time: ' num2str(OPTIONS(21)) ' seconds']); % GvW
disp(['f: ' num2str(f) ', max(g): ' num2str(mg)]); %GvW

% Auxilary function PSD, GvW
function [psd]=psd(nx)
if nx<0 || rem(nx,1)~=0;
  error('nx must be an integer > 0');
end;
psd=randn(nx); d=diag(abs(randn(1,nx))); psd=psd*d*psd';