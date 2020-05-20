  function [PHIk,GAMk,Qk,Rk,Mk,Vk,gk]=edortv(a,b,q,r,v,ts,z)
%EDORTV: Equivalent Discrete Optimal Control Problem Time-varying case
%
%    [PHIk,GAMk,Qk,Rk,Mk,Vk,gk]=edortv(a,b,q,r,v,ts,z)
%
%    function to calculate equivalent discrete optimal regulator matrices
%    for one sample interval [t ,t   ) for a linear time varying stochastic
%                              k  k+1
%    system A,B,W with quadratic cost J(Q,R) using the trapezium numerical
%    integration rule and a second order approximation of exp(A*dt)
%    where dt is the stepsize of the numerical integration.
%
%    Regulator problem:
%
%    x = A(t) x + B(t) u + W(t)
%    J = int(0,t )(x'Q(t)x+2*x'*Z(t)*u+u'R(t)u)+x(t )Hx(t )
%               f                                  f     f 
%
%    t , k=0,1,2,..,N-1: sampling instants, not necessarily equidistant
%     k
%
%    t = t
%     f   N
%
%    x(t ) = x :    Measured or reconstructed states
%       k     k
%
%    u(t)=u(t )=u , t <=t<t   , k=0,1,...,N-1:   Piecewise constant control
%            k   k   k     k+1
%
%    v = cov(W(t)):  Covariance matrix of white noise W
%    ni:  Number of numerical integration steps during [t ,t   )
%                                                        k  k+1
%    dt =(t   -t )/ni
%          k+1  k
%
%    Outputs:
%
%    PHIk,GAMk   : equivalent discrete system matrices
%    Qk,Rk,Mk  : equivalent discrete criterion matrices
%    Vk        : equivalent discrete noise covariance matrix
%    gk        : average cost due to white noise
%
%    Equivalent discrete regulator problem:
%
%    x    = PHI  x + GAM  u
%     k+1      k k      k  k
%
%    j = sum(k=0,N-1) (x'Q x  +  u' R u + 2*x' M u ) + x'Hx
%                       k k k     k  k k     k  k k     N  N
%
%    Inputs:
%
%    a         :[ A(t ) A(t +dt) A(t +2*dt) ... A(t   ) ]
%                    k     k        k              k+1
%    b         :[ B(t ) B(t +dt) B(t +2*dt) ... B(t   ) ]
%                    k     k        k              k+1
%    q         :[ Q(t ) Q(t +dt) Q(t +2*dt) ... Q(t   ) ]
%                    k     k        k              k+1
%    v         :[ V(t ) V(t +dt) V(t +2*dt) ... V(t   ) ]
%                    k     k        k              k+1
%    r         :[ R(t ) R(t +dt) R(t +2*dt) ... R(t   ) ]
%                    k     k        k              k+1
%
%         L.G. Van Willigenburg, W.L. De Koning, 10-01-02.
%
  if nargin==6; z=zeros(size(b));
  elseif nargin~=7;
    error('edortv requires 6 or 7 input arguments');
  end
  [nt,mt]=size(a);n=sqrt(nt);ni=mt-1;
  if rem(n,1)~=0;error('edortv: size mismatch argument 1');end;
  [nh,mh]=size(b);m=nh/n;mm=m*m;nm=n*m;
  if rem(m,1)~=0|mh~=mt;error('edortv: size mismatch argument 2');end;
  [nh,mh]=size(q);
  if nh~=nt|mh~=mt;error('edortv: size mismatch argument 3');end
  [nh,mh]=size(r);
  if nh~=mm|mh~=mt;error('edortv: size mismatch argument 4');end
  [nh,mh]=size(v);
  if nh~=nt|mh~=mt;error('edortv: size mismatch argument 5');end
  [nh,mh]=size(ts);
  if nh~=1|mh~=1;error('edortv: size mismatch argument 6');end
  [nh,mh]=size(z);
  if nh~=nm|mh~=mt;error('edortv: size mismatch argument 7');end

  Qk=zeros(n); Vk=zeros(n); GAMk=zeros(n,m); Mk=zeros(n,m);
  Rk=zeros(m); trud=0.0; gk=0.0;
  PHIk=eye(n); e=eye(n); delt=ts/ni; delta=0.5*delt;

  [bh]=getmat(b,n,1);[qh]=getmat(q,n,1); 
  [rh]=getmat(r,m,1); [vh]=getmat(v,n,1);
  [zh]=getmat(z,n,1);

  GAMkud=delta*bh; Qkud=delta*qh; Vkud=delta*vh; Rkud=delta*rh; Mkud=delta*zh;

  for i=1:ni 
    ind=i+1;
    [ah1]=getmat(a,n,i);[ah2]=getmat(a,n,ind);[bh]=getmat(b,n,ind);
    [qh]=getmat(q,n,ind);[vh]=getmat(v,n,ind);[rh]=getmat(r,m,ind);
    [zh]=getmat(z,n,ind);
    aud=delta*(ah1+ah2);
    PHIkud=e+aud+0.5*aud*aud;   PHIk=PHIkud*PHIk;
    GAMk=PHIkud*(GAMk+GAMkud);   GAMkud=delta*bh;                 GAMk=GAMk+GAMkud;
    Vk=PHIkud*(Vk+Vkud)*PHIkud'; Vkud=delta*vh;                   Vk=Vk+Vkud;
    Qk=Qk+Qkud;                  Qkud=delta*PHIk'*qh*PHIk;        Qk=Qk+Qkud;
    gk=gk+trud;                  trud=delta*trace(Vk*qh);         gk=gk+trud;
    Mk=Mk+Mkud;                  Mkud=delta*PHIk'*(qh*GAMk+zh);   Mk=Mk+Mkud;
    h=GAMk'*zh;                  h=h+h';
    Rk=Rk+Rkud;                  Rkud=delta*(rh+GAMk'*qh*GAMk+h); Rk=Rk+Rkud;
  end