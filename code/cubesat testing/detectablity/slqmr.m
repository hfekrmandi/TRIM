% SLQMR: Function computes [Q M; M' R] or F of LQG problem based on
%        scaled second order terms (second derivatives)
%
% 1      [Q M; M' R] based on second deriv. of f(x,u,t), g(x,u,t) and possibly H(x,u,t)
%        The latter requires la, the costate, to be supplied as an input
% 
% 2      F based on second deriv. of phi(x(tf))
%
%        [qmr,hest,qmrs,hests]=slqmr(fname,dims,t,x,u,dxu,ds,str,la,varargin);
%
%        Input:
%        fname : name f,g,phi
%        dims  : dimensions
%        t     : time
%        x     : state value
%        u     : control
%        ds    : vector with scaling to 1 factors of x,u,y
%        str   : string either 'f','g' or 'phi'
%        la    : co-state (optional used to include H(x,u,t))
%        dxu   : vector with perturbations for x,u
%                dx(k)=dxu(k), k=1:length(x)
%                du(k)=dxu(length(x)+k), k=1:length(u)
%   varargin   : additional inputs. varargin{1:2}={[pv] [extinp]}
%
%        Output:
%        qmr   : [Q M;M' R] of F
%        hest  : Hessian(k), k=1:nx
%        qmrs  : Scaled to 1 diag([Q M;M' R])
%        hests : Scaled to 1 Hessian(k), k=1:nx
%        

function [qmr,hest,qmrs,hests]=slqmr(fname,dims,t,x,u,dxu,ds,str,la,varargin)
                         
dimx=dims(1); dimu=dims(2); dimy=dims(3); dimd=dims(4); dimp=dims(5); dimxd=dimx+dimd;

% Check dimensions of inputs
[nx,mx]=size(x); [nu,mu]=size(u); [nxu]=size(dxu,1); [nd,md]=size(ds);
if nxu~=nx+nu; error(' Dimensions x,u,dxu incompatible'); end;
if nd~=dimxd+dimu+dimy+1; error(' Dimension ds incompatible'); end;

[nh,mh]=size(ds);
if mh==nh; ds=diag(ds); end;
if min(ds)<=0; error(' ds(i) > 0 violated'); end;

if nargin<6; error(' At least 6 input arguments required'); end
if isempty(u); u=zeros(0,1); end;
if isempty(la); la=zeros(nx,1); end;

% Determine state and control deviations for derivative computations
if nx~=dimxd+1 && nx~=dimxd
  error(' nx incompatible')
end
  
% Obtain mt (method) to compute QMR
[wf,meth,mt]=lqgwfgi(str);

% Scale lambda
if str=='f'
  if nx==dimx; la(1:dimx)=la(1:dimx)./ds(1:dimx);
  else la(1:dimx+1)=la(1:dimx+1)./ds(1:dimx+1);
  end
end  

flag=1; %keyboard

% Compute f(x,u)
[xdot]=feval(fname,t,x,u,varargin{:}); [nf,mf]=size(xdot);

% Compute Hessians, Scaled Hessians & QMR
% Initialize
qt=[]; qts=[]; hest=[]; hests=[]; qmrs=zeros(nxu); qmrsh=qmrs;

% Use Automatic differentiation if intlab is installed
if exist('hessianinit','file')
% Computation by automatic differentiation using intlab
  xuH=hessianinit([x;u]);
  xdH=feval(fname,t,xuH(1:nx),xuH(nx+1:nxu),varargin{:});
  sMat=1./ds(1:nxu); sMat=sMat*sMat';

% Collect Hessians, scaled Hessians & compute QMR
  for k=1:nf
  % Obtain and store unscaled Hessians
    q=full(xdH(k).hx); hest=[hest q];

  % Compute and store scaled Hessians
    qs=ds(k)*sMat.*q; hests=[hests qs];

  % Compute symmetric nonnegative approximation qsnn of qs
    [u,d,v]=svd(qs); qsnn=u*d*u';
    [qsnn]=cin2nn(qs,mt);

  % Obtain scaling factors and method to compute QMR
    [wfgp,meth]=lqgwfgi(str); wp=lqgwfgi('phi');

  % Build scaled to 1 QMR two methods possible
    if str=='f'; nfgp=nf-1; else nfgp=nf; end;
    if meth==0
      if str=='f'
        if k==nf; w=wp*max(abs(la(k)),1);
        else w=max(wp*abs(la(k)),wfgp/nfgp);
        end
      else w=wfgp/nfgp;
      end;
      qmrs=qmrs+w*qsnn;
    else
      if str=='f'
        if k==nf; w=wp; else w=wfgp/nfgp; end;
        qmrsh=qmrsh+w*la(k)*qs;
      else
        w=wfgp/nfgp; qmrs=qmrs+w*qsnn;
      end;
    end
  end
  
elseif exist('fmad','file')

  % Computation by automatic differentiation using tomlab
  xuH=fmad([x;u],eye(nxu)); xuH=fmad(xuH,eye(nxu));
  xdH=feval(fname,t,xuH(1:nx),xuH(nx+1:nxu),varargin{:});
  sMat=1./ds(1:nxu); sMat=sMat*sMat';

% Collect Hessians, scaled Hessians & compute QMR
  for k=1:nf
  % Obtain and store unscaled Hessians
    q=getinternalderivs(getinternalderivs(xdH(k)));
    if isempty(q); q=zeros(nxu); end; hest=[hest q];
  % Compute and store scaled Hessians
    qs=ds(k)*sMat.*q; hests=[hests qs];

  % Compute symmetric nonnegative approximation qsnn of qs
  %  [u,d,v]=svd(qs); qsnn=u*d*u';
    [qsnn]=cin2nn(qs,mt);

  % Obtain scaling factors and method to compute QMR
    [wfgp,meth]=lqgwfgi(str); wp=lqgwfgi('phi');

  % Build scaled to 1 QMR two methods possible
    if str=='f'; nfgp=nf-1; else nfgp=nf; end;
    if meth==0
      if str=='f'
        if k==nf; w=wp*max(abs(la(k)),1);
        else w=max(wp*abs(la(k)),wfgp/nfgp);
        end
      else w=wfgp/nfgp;
      end;
      qmrs=qmrs+w*qsnn;
    else
      if str=='f'
        if k==nf; w=wp; else w=wfgp/nfgp; end;
        qmrsh=qmrsh+w*la(k)*qs;
      else
        w=wfgp/nfgp; qmrs=qmrs+w*qsnn;
      end;
    end
  end
  
else

  % Computation by finite differences if intlab not installed
  for i=1:nxu
    for j=i:nxu

    % Compute df/dxui, dxuj>0
      dx=zeros(nx,1); du=zeros(nu,1);
      if j<=nx; dx(j)=dxu(j); else du(j-nx)=dxu(j); end;

      if i<=nx; dx(i)=dxu(i); else du(i-nx)=dxu(i); end;
      [xdotf]=feval(fname,t,x+dx,u+du,varargin{:});
      if i<=nx; dx(i)=-dxu(i); else du(i-nx)=-dxu(i); end;
      [xdotb]=feval(fname,t,x+dx,u+du,varargin{:});

      if i~=j

      % Store df/dxui, dxuj>0
        dff=(xdotf-xdotb)/(2*dxu(i));

      % Compute df/dxui, dxuj<0
        if j<=nx; dx(j)=-dxu(j); else du(j-nx)=-dxu(j); end;

        if i<=nx; dx(i)=dxu(i); else du(i-nx)=dxu(i); end;
        [xdotf]=feval(fname,t,x+dx,u+du,varargin{:});
        if i<=nx; dx(i)=-dxu(i); else du(i-nx)=-dxu(i); end;
        [xdotb]=feval(fname,t,x+dx,u+du,varargin{:});

      % Store df/dxui, dxuj<0
        dfb=(xdotf-xdotb)/(2*dxu(i));

      % Compute d2f/dxuidxuj
        h=(dff-dfb)/(2*dxu(j));
      else
      % Compute d2f/dxuidxuj, add it to qt, scale h add it to qts
        h=(xdotf-2*xdot+xdotb)/(dxu(i)*dxu(i));
      end

  %   Add h to qt
      qt=[qt h]; 
  %		Scale h and add to qts
      h=h/(ds(i)*ds(j));
      if str=='f'; h=ds(1:nf).*h;
      elseif str=='g';
      % Output indices and scaling of h
        yind=dimxd+dimu+2:dimxd+dimu+dimy+1;
        h=ds(yind).*h;
      elseif str=='phi'; h=ds(dimx+1)*h;
      else
        error('  str should be ''f'',''g'' or ''phi''');
      end
      qts=[qts h]; %if max(max(abs(h)))~=0; keyboard; end
    end
  end

  % qt(k,:) contains d2f(k)/dxuidxuj, k=1:nx
  % the lower diagonal part of q(k)

  % Create matrix index array
  ind=[];
  for k=1:nxu
    ind=[ind (k-1)*nxu+k:k*nxu];
  end

  % Build qmr and hest
  for k=1:nf

  % Clear and create qh = q(:)
    qh=zeros(nxu*nxu,1); qhs=qh;
    qh(ind)=qt(k,:)'; qhs(ind)=qts(k,:)';

  % Clear and create lower diagonal
  % matrix q, q(:) = qh 
    q=zeros(nxu); qs=q;
    q(:)=qh; qs(:)=qhs;

  % Symmetrize q and add q to hessians
    q=(q+q')-diag(diag(q));  qs=(qs+qs')-diag(diag(qs));
    hest=[hest q]; hests=[hests qs];

  % Compute symmetric nonnegative approximation qsnn of qs
  %  [u,d,v]=svd(qs); qsnn=u*d*u';
    [qsnn]=cin2nn(qs,mt);

  % Obtain scaling factors and method to compute QMR
    [wfgp,meth]=lqgwfgi(str); wp=lqgwfgi('phi');

  % Build scaled to 1 QMR two methods possible
    if str=='f'; nfgp=nf-1; else nfgp=nf; end;
    if meth==0
      if str=='f'
        if k==nf; w=wp*max(abs(la(k)),1);
        else w=max(wp*abs(la(k)),wfgp/nfgp);
        end
      else w=wfgp/nfgp;
      end;
      qmrs=qmrs+w*qsnn;
    else
      if str=='f'
        if k==nf; w=wp; else w=wfgp/nfgp; end;
        qmrsh=qmrsh+w*la(k)*qs;
      else
        w=wfgp/nfgp; qmrs=qmrs+w*qsnn;
      end;
    end
  end
end
% End AD/FD computation

% Addition needed for second method
if str=='f' & meth~=0; qmrs=qmrs+cin2nn(qmrsh,mt); end;

% Scale back qmr
qmr=diag(ds(1:nxu))*qmrs*diag(ds(1:nxu));