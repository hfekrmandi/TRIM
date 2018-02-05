% DPRLUDUB: One iteration of two Lyapunov equations that determine the
%           deterministic parameter optimal reduced order LQG compensator
%           based on UDU factorizations.
%       
%           L.G. van Willigenburg, W.L. de Koning 13-1-'99.
%
%           Input:  p,g,c,q,r,v,w; LQG problem parameters.
%                   pp; Augmented system second moment matrix.
%                   ss; Lagrange multipliers.
%
%           Output: pp; Augmented system second moment matrix.
%                   ss; Lagrange multipliers.
%                   f,k,l: compensator matrices.
%                   sigp, sigs: compensator costs.
%                   tps       : trace(pt+st).
%
%           See also dprotil
%

  function [pp,ss,f,k,l,sigp,sigs,tps,pt,st]=dprludub(pp,ss,p,g,c,q,r,v,w,mc,me,b)

  [nxn,nx]=size(p); [nt,nt]=size(pp);
  [nu,nu]=size(r); [ny,ny]=size(w);
  nc=nt-nx; nx1=nx+1;

  if exist('eps')==1;
    if eps>0 & eps<1e-6; epsinv=0.1*max(nx,nxn)*eps;
    else; epsinv=1e-12; end
  else; epsinv=1e-12; end

  if nt==nx;
    f=[]; l=[]; pf=1;
  else
    pf=0;
  end;

  p1=pp(1:nx,1:nx); p12=pp(1:nx,nx1:nt);
  p2=pp(nx1:nt,nx1:nt); pt=ud2sym(p1);

  if min(size(r))>0;
    ru=sym2ud(r,epsinv,-epsinv);
  end
  if min(size(w))>0;
    wu=sym2ud(w,epsinv,-epsinv);
  end

  [nt,nt]=size(ss); ncn=nt-nxn; nxn1=nxn+1;

  if nt==nxn
    f=[]; k=[]; sf=1;
  else
    sf=0;
  end

  s1=ss(1:nxn,1:nxn); s12=ss(1:nxn,nxn1:nt);
  s2=ss(nxn1:nt,nxn1:nt); st=ud2sym(s1);

  if pf==1

    [k0]=ktudu(p1,wu,c,[],epsinv,-epsinv);
    [k0]=udinv(k0,epsinv,-1e12); k0=(me+p*pt*c')*udt2sym(k0);

    if ncn~=0
      [uh]=udinv(s2,epsinv,-1e12); [uh]=ud2ud(uh); hh=-uh'*s12'; k=hh*k0;
      kwk=ud2sym(ktudu(wu,zeros(ncn),k,[],epsinv,-epsinv));   
      pa=[p; k*c]; qa=q; h=me*k'; va=[v  h;  h' kwk];
    else
      pa=p; qa=q; va=v;
    end

  elseif sf==1

    [l0]=ktudu(s1,ru,g',[],epsinv,-epsinv);
    [l0]=udinv(l0,epsinv,-1e12); l0=udt2sym(l0)*(g'*st*p+mc');

    if nc~=0
      [uh]=udinv(p2,epsinv,-1e12); [uh]=ud2ud(uh); gg= uh'*p12'; l=l0*gg';
      lrl=ud2sym(ktudu(ru,zeros(nc),l',[],epsinv,-epsinv));
      pa=[p -g*l]; va=v; h=-mc*l; qa=[q h ; h' lrl];
    else
      pa=p; va=v; qa=q;
    end

  elseif min(size(r))==0

    [k0]=ktudu(p1,wu,c,[],epsinv,-epsinv);
    [k0]=udinv(k0,epsinv,-1e12); k0=(me+p*pt*c')*udt2sym(k0);

    [uh]=udinv(s2,epsinv,-1e12); [uh]=ud2ud(uh); hh=-uh'*s12';
    [uh]=udinv(p2,epsinv,-1e12); [uh]=ud2ud(uh); gg= uh'*p12';

    k=hh*k0; f=hh*(p-k0*c)*gg'; l=[];

    kwk=ud2sym(ktudu(wu,zeros(ncn),k,[],epsinv,-epsinv));

    pa=[p zeros(nxn,nc); k*c f];
    qa=[q zeros(nx,nc); zeros(nc,nx+nc)]; 
    h=me*k'; va=[v  h;  h' kwk];

  elseif min(size(w))==0

    [l0]=ktudu(s1,ru,g',[],epsinv,-epsinv);
    [l0]=udinv(l0,epsinv,-1e12); l0=udt2sym(l0)*(g'*st*p+mc');

    [uh]=udinv(s2,epsinv,-1e12); [uh]=ud2ud(uh); hh=-uh'*s12';
    [uh]=udinv(p2,epsinv,-1e12); [uh]=ud2ud(uh); gg= uh'*p12';

    l=l0*gg'; f=hh*(p-g*l0)*gg'; k=[];

    lrl=ud2sym(ktudu(ru,zeros(nc),l',[],epsinv,-epsinv));

    pa=[p -g*l; zeros(ncn,nx) f];
    h=-mc*l; qa=[q h ; h' lrl];
    va=[v zeros(nxn,ncn); zeros(ncn,nxn+ncn)];

  else

    [l0]=ktudu(s1,ru,g',[],epsinv,-epsinv);
    [l0]=udinv(l0,epsinv,-1e12); l0=udt2sym(l0)*(g'*st*p+mc');
    [k0]=ktudu(p1,wu,c,[],epsinv,-epsinv);
    [k0]=udinv(k0,epsinv,-1e12); k0=(me+p*pt*c')*udt2sym(k0);

    [uh]=udinv(s2,epsinv,-1e12); [uh]=ud2ud(uh); hh=-uh'*s12';
    [uh]=udinv(p2,epsinv,-1e12); [uh]=ud2ud(uh); gg= uh'*p12';

    f=hh*(p-k0*c-g*l0)*gg'; k=hh*k0; l=l0*gg';

    lrl=ud2sym(ktudu(ru,zeros(nc),l',[],epsinv,-epsinv));
    kwk=ud2sym(ktudu(wu,zeros(ncn),k ,[],epsinv,-epsinv));

    pa=[p -g*l; k*c f];
    h=-mc*l; qa=[q h ; h' lrl];
    h=me*k'; va=[v h ; h' kwk];
  end

  tps=trace(pt)+trace(st);
  sigp=trace(qa*ud2sym(pp)); sigs=trace(va*ud2sym(ss));
  qa=sym2ud(qa,epsinv,-1e12); va=sym2ud(va,epsinv,-1e12);

%  pp=ktudu(pp,va,pa ,[],epsinv,-1e12);
  ss=ktudu(ss,qa,pa',[],epsinv,-1e12);
