% UD2SYM : Recover P from UD an UDU factorization of P
%
%          See Bierman 1977: Factorization Methods for
%          discrete sequential estimation
%
%          [p]=ud2sym(ud)
%          p=u*d*u'; u unit upper triangular, d diagonal
%                    both stored in ud
%
%          GvW 7-1-'99
%
  function [p]=ud2sym(ud)
  
  [n,m]=size(ud);
  if n~=m; error('ud must be square'); end;
  if n==0; error('  Compatible but empty inputs'); end;

  d=diag(diag(ud));
  u=ud-d+diag(ones(1,n));
  p=u*d*u';

