% PS2PHSH: Convert pp,ss (Lyapunov equations)
%          to p,s,ph,sh (Projection equations)
%
%          [pst,psht]=ps2phsh(ppsst,nx,nc);
%
%      Input:
%          ppsst  : [ps(1)...ps(N)]
%          nx   : dimension system
%          nc   : dimension compensator
%
%      Output:
%          pst  : [ps(1)....ps(N)]
%          psht : [psh(1)...psh(N)]
%
%          L.G. Van Willigenburg, W.L. De Koning 27-06-'96
%
  function [pst,psht]=ps2psh(ppsst,nx,nc)

  pst=[]; psht=[];
  [n,m]=size(ppsst);

  if min(n,m)>0;

    nxc=nx+nc; nx1=nx+1;
    if n~=nxc; error(' Incompatible dimensions input'); end;
    N=round(m/nxc);
  
    for i=1:N
      pn2=i*nxc; pn1=pn2-nxc+1; ppss=ppsst(:,pn1:pn2);
      if nc==0; ps=ppss; psh=zeros(nx);
      else
        ps1=ppss(1:nx,1:nx); ps12=ppss(1:nx,nx1:nxc);
        ps2=ppss(nx1:nxc,nx1:nxc);
        % pips2=pinv(ps2); psh=ps12*pips2*ps12'; ps=ps1-psh; %old
        % Symmetrize, svd, pseudo inverse % new
        ps1=0.5*(ps1+ps1'); ps2=0.5*(ps2+ps2');
        [u,d,v]=svd(ps2); d=diag(d); ld=length(d); dp=sum(d>ld*eps);
        ds=[1./sqrt(d(1:dp)); zeros(ld-dp,1)];
        pips2=u*diag(ds); psh=ps12*pips2; psh=psh*psh'; ps=ps1-psh;
      end
      pst=[pst ps]; psht=[psht psh];
    end
  end
