% UDT2SYM: Recover P=U'DU, UD an UDU factorization of P
%
%          See Bierman 1977: Factorization Methods for
%          discrete sequential estimation
%
%          [p]=udt2sym(ud)
%          p=u'*d*u, used to compute inv(p):
%          inv(p)=udt2sym(udinv(sym2ud(p)))
%
%          see also udinv, ud2sym, sym2ud
%          GvW 7-1-'99
%
  function [p]=udt2sym(ud)
  
  [n,m]=size(ud);
  if n~=m; error('ud must be square'); end;
  if n==0; error('  Compatible but empty inputs'); end;

  d=diag(diag(ud));
  u=ud-d+diag(ones(1,n));
  p=u'*d*u;

