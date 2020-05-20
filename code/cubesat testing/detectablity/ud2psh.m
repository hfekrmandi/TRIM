% UD2PSH : Convert ud factorisation of PP, SS to P,Phat or S, Shat
%
%  function [ps,psh]=ud2psh(udps,nc);
%
  function [ps,psh]=ud2psh(udps,nc);

  [nx1,nx2]=size(udps);

  if nx1~=nx2; error(' udps must be square upper triangular'); end;
  if nc>=nx1; error(' nc must be smaller than size udps'); end;

  nx=nx1-nc; nx1=nx+1;

  ps=ud2sym(udps(1:nx,1:nx));

  if nc~=0
    u12=udps(1:nx,nx1:nx2); u2=udps(nx1:nx2, nx1:nx2);
    [u2,d2]=ud2ud(u2);
    psh=u12*d2*u12';
  else
    psh=zeros(nx);
  end;
