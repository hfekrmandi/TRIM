% GENNOI : Generate measurement, system and initial state noise 
%
%          [x0,dx0,vh,wh]=gennoi(dims,lqgarr,tk,x0,wt,cw,mulnoi);
%
%          Input:
%
%             dims    : dimensions
%             lqgarr  : array with equivalent discrete LQG problem data
%             tk      : sampling instants
%             dt      : time step numerical integration
%             x0      : initial state non-linear system
%             vt      : continuous-time white noise intensity matrices
%             wt      : discrete-time white noise covariance matrices
%             cw~=0   : continuous white system noise generation else discrete-time
%             mulnoi  : Noise multiplication factor
%
%          Output:
%
%             x0      : Extended non-linear system state
%             dx0     : Realization initial state error
%             vh      : Relaization system noise, either continuous or discrete-time
%             wh      : Realization discrete-time measurement noise
%
% GvW 16-5-2005

  function [x0,dx0,vh,wh]=gennoi(dims,lqgarr,tk,dt,x0,vt,wt,cw,mulnoi)
  
  dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;
  
  [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me,x0m,x0c]=getedo(lqgarr,1);
  
% Specify the initial compensator state errors
  [u,d,u]=svd(x0c); x0c=u*diag(sqrt(diag(d)));
  dx0=mulnoi*[x0c*randn(length(x0)-1,1)]; dx0=adj2rc(dx0,dimxd);

% Measurement (wh) and sytem noise (vh) generation
  ic=0; N=lqgarr(1);
  wh=zeros(ldi,N); vh=[];
  for i=1:N-1
    [nx,nxn,mi,ldi,p,g,q,r,mc,v,c,w,me]=getedo(lqgarr,i);
    if ~cw
    % Discrete-time white system noise generation vh
      [u,d,u]=svd(v); vs=u*diag(sqrt(diag(d)));
      vs=mulnoi*vs*randn(nx,1); vs=adj2rc(vs,dimxd,dimd);
      vh=[vh vs];
    else
    % Continuous-time white system noise generation vh
      j=round((tk(i+1)-tk(i))/dt);
      for k=1:j; 
        vs=zeros(dimxd); vs(:)=vt(:,ic+k);
        [u,d,u]=svd(vs); vs=u*diag(sqrt(diag(d)));
        vs=mulnoi*vs*randn(nx); vs=adj2rc(vs,dimxd);
        vh=[vh vs];
      end;
      ic=ic+j;
    end
    % Discrete-time white measurement noise generation wh
    ws=zeros(ldi); ws(:)=wt(:,i); [u,d,u]=svd(ws);
    ws=u*diag(sqrt(diag(d))); wh(:,i)=mulnoi*ws*randn(ldi,1);
  end