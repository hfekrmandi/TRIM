% UDUADD : Add (scaled) UDU factored matrices P,V: P=pc*P+vc*V, pc,vc scalars
%
%          See Bierman 1977: Factorization Methods for
%          discrete sequential estimation pp. 132-133
%
%          [udp,r]=uduadd(udp,udv,pc,vc,tolp,toln)
%
%          input:
%          udp : ud of p
%          udv : ud of v
%          pc  : scale factor P, default 1
%          vc  : scale factor V, default 1
%
%          output:
%          udp: ud of pn
%          r  : rank of d
%
  function [udp,r]=uduadd(udp,udv,pc,vc,tolp,toln)

  if nargin<2; error(' 2-6 input arguments'); end;
  if nargin<3; pc=1; vc=1; tolp=1e-12; toln=tolp;
  elseif nargin<4
    [n,m]=size(pc);
    if n>1 | m>1; error(' pc must empty or scalar'); end;
    if n==0 | m==0; pc=1; end;
    vc=1; tolp=1e-12; toln=tolp;
  elseif nargin<5;
    [n,m]=size(vc);
    if n>1 | m>1; error(' vc must empty or scalar'); end;
    if n==0 | m==0; vc=1; end;
    tolp=1e-12; toln=tolp;
  elseif nargin<6;
    toln=tolp;
  elseif nargin>6;
    error(' 2-6 input arguments');
  end;

  tolp=max(0,tolp);
  if (toln>tolp); error('  toln > tolp'); end;

  [n,m]=size(udp); [n1,m1]=size(udv);
  if n~=m | n1~=n | m1~=m; error('  udp, udv must be square of equal size'); end;
  if n==0 | m==0; error('  Compatible but empty inputs'); end;

  nud=n+n; ud=[udp udv];

  a=zeros(1,nud); v=zeros(1,nud); d=zeros(1,nud);

  for i=1:n; d(i)=pc*ud(i,i); ud(i,i)=1; end;

  for i=1:nud-n; d(i+n)=vc*ud(i,n+i); ud(i,n+i)=1; end;

  r=n;
  for j=n:-1:1
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

  for j=2:n;
    for i=1:j-1;
      ud(i,j)=ud(j,i); ud(j,i)=0;
    end;
  end;
  udp=ud(1:n,1:n);
