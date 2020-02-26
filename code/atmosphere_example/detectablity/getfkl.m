% GETFKL : Get or display compensator matrices f,k,l
%          for time i or all times i=1,2,..,N
%          from the compensator array xfkl
%          Time shift i=0,1,..,N -> i=1,2,..,N+1
%
%          [nc,ncn,nu,ny,f,k,l,x0]=getfkl(xfkl,i,tol);
%
%	   input:
%          xfkl   : compensator in array format
%	   Array format:
%	         [N;ind(1);..;ind(N);
%	         nc(1);nc(2);nu(1);ny(1);x0(:);
%	         f(1)(:);k(1)(:);l(1)(:)
%	         ...
%	         nc(i);nc(i+1);nu(i);ny(i);
%	         f(i)(:);k(i)(:);l(i)(:)];
%	         ...
%	         nc(N-1);nc(N);nu(N-1);ny(N-1);
%	         f(N-1)(:);k(N-1)(:);l(N-1)(:)];
%          i      : Time index 1,2,..,N
%                   if not specified the full compensator is displayed
%
%    output :
%          If no output specified display occurs
%          nc,ncn,nu,ny : nc(i),nc(i+1),nu(i),ny(i)
%          f,k,l        : f(i), k(i), l(i)
%          x0           : x(1)
% 
  function [nc,ncn,nu,ny,f,k,l,x0]=getfkl(xfkl,i,tol);

  if nargin==1; i=0; tol=1e-60;
  elseif nargin==2; tol=1e-60;
  elseif nargin>3;
    error(' One to three input arguments required');
  end
  
  if nargout~=0 & nargin==1; error(' Time i must be specified'); end;
  [nh,mh]=size(xfkl);
  if mh~=1; error(' First input argument must be a column vector'); end;

% Determine N (=N+1) and check i
  N=xfkl(1); oi=1;

  if i==0; dall=1; ih=1:N-1;
  elseif isempty(i); dall=1; ih=1:N-1;
  else; 
    dall=0; ih=i; 
    if i>N-1 | i<1; error(' 1 <= i <= N-1 violated'); end;
  end

% Obtain x0
  oi=xfkl(2); n=xfkl(oi+1); x0=zeros(n,1);
  x0(:)=xfkl(oi+4+1:oi+4+n);
  if nargout==0; disp('x(0)'''); disp(roundmat(x0',tol)); end;

  for i=ih
  % Determine the offset o(i)
    oi=xfkl(1+i);

  % Determine n(i),n(i+1),m(i),l(i) and
    oi=oi+1; nc=xfkl(oi); oi=oi+1; ncn=xfkl(oi);
    oi=oi+1; nu=xfkl(oi); oi=oi+1; ny=xfkl(oi);

  % Extra offset to skip x0 if i==1
    if i==1; oi=oi+nc; end;

% Obtain matrices
    n=ncn;  m=nc; nm=n*m;
    f=zeros(n,m); f(:)=xfkl(oi+1:oi+nm); oi=oi+nm;
    n=ncn; m=ny; nm=n*m;
    k=zeros(n,m); k(:)=xfkl(oi+1:oi+nm); oi=oi+nm;
    n=nu; m=nc; nm=n*m;
    l=zeros(n,m); l(:)=xfkl(oi+1:oi+nm); oi=oi+nm;
    if xfkl(oi+1)~=ncn; error(' ncn check error'); end;

    if nargout==0
      is=int2str(i-1); disp(' ');
      disp(['  F(' is ')']); disp(roundmat(f,tol));
      disp(['  K(' is ')']); disp(roundmat(k,tol));
      disp(['  L(' is ')']); disp(roundmat(l,tol));
    end
  end