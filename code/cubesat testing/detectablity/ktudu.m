% KTUDU :  Kalman time update based on UDU factorizations
%          of symmetric p and v
%
%          See Bierman 1977: Factorization Methods for
%          discrete sequential estimation pp. 132-133
%
%          time update:
%          pn=pd*p*pd'+v, xn=pd*x;
%
%          [udp,r,x]=ktudu(udp,udv,pd,x,tolp,toln)
%
%          input:
%          udp : ud of p
%          udv : ud of v
%          pd  : may not be square but then x will be empty
%          x   : state estimate
%
%          output:
%          udp: ud of pn
%          x  : time updated state estimate
%          r  : rank of d
%
  function [udp,r,x]=ktudu(udp,udv,pd,x,tolp,toln)

  if nargin<3; error(' 3-6 input arguments'); end;
  if nargin<=4; tolp=1e-12; toln=tolp;
  elseif nargin==5; toln=tolp; end; tolp=max(0,tolp); 
  if (toln>tolp); error('  toln > tolp'); end;

  [n1,n2]=size(udp); [m1,m2]=size(udv); [m,n]=size(pd);
  if n1~=n | n2~=n; error('  udp must be square and compatible with pd'); end;
  if m1~=m | m2~=m; error('  udv must be square and compatible with pd'); end;

  if nargin>3; [nx,mx]=size(x);
    if nx~=0 & mx~=0; xflag=1; else; xflag=0; end;
  else; xflag=0; end;

  if xflag==1
    if m~=nx;  error(' size pd and x do not match'); end;
    if n~=m; error('  udp, pd, and x should have equal n'); end;
  end;

  if n==0 | m==0; error(' Compatible but empty inputs'); end;

  nud=n+m;

  if n>m; ud=[udp [udv; zeros(n-m,m)]];
  else; ud=[[udp; zeros(m-n,n)] udv];
  end;

  a=zeros(1,nud); v=zeros(1,nud); d=zeros(1,nud);
  if xflag==1; v(1)=pd(1,1)*x(1); end;

  for j=n:-1:2
    for i=1:j; d(i)=ud(i,j); end %5
    for i=1:m
      ud(i,j)=pd(i,j);
      for k=1:j-1
        ud(i,j)=ud(i,j)+pd(i,k)*d(k);
%        disp([i j k]); disp([ud(i,j) pd(i,k) d(k)]);
      end
    end %10

    if xflag==1;
      v(j)=0; for k=1:n; v(j)=v(j)+pd(j,k)*x(k); end %15
      v(1)=v(1)+pd(1,j)*x(j);
    end %if
  end; %20
  d(1)=ud(1,1);
      
% v=pd*x, u=pd*u and diagonals retrieved

  for j=1:m
    if xflag==1; x(j)=v(j); end; ud(j,1)=pd(j,1);
  end %30

% x=pd*x, ud=pd*ud completed

%  keyboard;
  for i=1:nud-n; d(i+n)=ud(i,n+i); ud(i,n+i)=1; end;

  r=n;
  for j=m:-1:1
    sig=0;
    for k=1:nud
      v(k)=ud(j,k); a(k)=d(k)*v(k); sig=sig+v(k)*a(k);
    end %40
    ud(j,j)=sig;
    if sig<toln; error(' toln violated');
    elseif sig<tolp; sig=0; dinv=0; r=r-1;
    else; dinv=1/sig; end;
    if j~=1
      jm1=j-1;
      for k=1:jm1
        sig=0;
        for i=1:nud
          sig=sig+ud(k,i)*a(i);
        end %50
        sig=sig*dinv;
        for i=1:nud
          ud(k,i)=ud(k,i)-sig*v(i);
          ud(j,k)=sig;
        end %60
      end %70
    end %if
  end %100

  for j=2:m;
    for i=1:j-1;
      ud(i,j)=ud(j,i); ud(j,i)=0;
    end;
  end;
  udp=ud(1:m,1:m);
