  function diglqg(lqgmatr,edorp,digcomp,xfkldat,compv4,mro,opt)
  % DIGLQG : Function that computes a digital LQG compensator
  %  
  %          diglqg(lqgmatr,edorp,digcomp,xfkldat,compv4)
  %
  %          Input:
  %          lqgmatr : Name of file containing the LQG problem data
  %          edorp   : Name of file containing the equivalent discrete-time LQG
  %                    problem data
  %
  %          Input/Output:
  %          digcomp : Name of file containing the digital LQG compensator data
  %          xfkldat : Name of ascii file with the digital LQG compensator data
  %          compv4  : Name of ascii file for optimal LQG compensation using Matlab 5
  %          mro     : Method to compute reduced-order LQG compensator
  %                    1 : Optimal Projection equations (defaul_t), function dprotva
  %                    2 : Lyapunov equations, function dprotvla
  %                    3 : UDU factored Lyapunov equations, function dprtvlau
  %          opt     : options for
  %                    1 : dprotva
  %                    2 : dprotvla
  %                    3 : dprtvlau
  %
  % GvW 18-7-2003

% Check inputs
  if nargin<5; error(' At least 5 inputs required'); end
  if ~isstr(lqgmatr); error(' 1st input must be a string'); end
  if ~isstr(edorp); error(' 2nd input must be a string'); end
  if ~isstr(digcomp); error(' 3rd input must be a string'); end
  if ~isstr(xfkldat); error(' 4th input must be a string'); end
  if ~isstr(compv4); error(' 5th input must be a string'); end
   
  if ~exist('mro','var'); mro=1; end;
  
% load equivalent discrete-time LQG problem data
  load(edorp); ns=length(tk)-1;

% Read prescribed compensator dimensions
  [nc]=lqgnc(lqgmatr,dims,ns,zeros(size(x0)));

% Compute the digital LQG compensator
  if ~exist('opt','var') || isempty(opt);
    opt=[]; opt(4)=10; opt(6)=1; opt(1)=1e-6; %opt(5)=0.4;
  end
  
  %av=[0 0.2 0.4 0.5 0.55 0.6 0.625:0.025:1]; opt=[opt av];
  if mro==1; [xfkl,sigp,sigs,ptt,stt,pht,sht]=dprotva(lqgarr,nc,opt);
  elseif mro==2 [xfkl,sigp,sigs,ppt,sst]=dprotvla(lqgarr,nc,opt);
  else [xfkl,sigp,sigs,ppt,sst,nsweeps]=dprtvlau(lqgarr,nc,opt); end
  sigp=sigp+costed,sigs=sigs+costed
  %sigp,sigs,costed

  if mro==1
    eval(['save ' digcomp ' sysdyn sysout phi xfkl ptt stt sigp sigs lqgarr tt qt rt mct lat vt wt']);
  else
    eval(['save ' digcomp ' sysdyn sysout phi xfkl ppt sst sigp sigs lqgarr tt qt rt mct lat vt wt']);
  end
  eval(['save -ascii ' xfkldat ' xfkl']);
  
  % Save data in file for control helicopter
  ncm=max(nc); f_t=[]; k_t=[]; l_t=[]; tol=1e-6;
  for i=1:ns
    [nci,ncn,nu,ny,fi,ki,li,x0]=getfkl(xfkl,i,tol);
    fc=zeros(ncm); fc(1:ncn,1:nci)=fi;
    kc=zeros(ncm,ny); kc(1:ncn,:)=ki;
    lc=zeros(nu,ncm); lc(:,1:nci)=li;
    f_t=[f_t fc]; k_t=[k_t kc]; l_t=[l_t lc];
  end
  
  n=nc; m=nu; l=ny; ts=satimes(1,2)-satimes(1,1);
  uhstar=ustar; ustar=uhstar(1:end-1,:);
  eval(['save ' compv4 ' n m l ts ustar ystar f_t k_t l_t -v4']);
%  save m:\comp n m l ts ustar ystar f_t k_t l_t -v4
  
  disp('Digital LQG compensator stored');