% UD2UD : Recover separate u,d from ud, an udu factorization
%         [u,d]=ud2ud(ud)
%
%         GvW 7-1-'99
%
  function [u,d]=ud2ud(ud)

  [n,m]=size(ud);
  if n~=m; error('ud must be square'); end;
  if n==0; error('  Compatible but empty inputs'); end;

  d=diag(diag(ud));
  u=ud-d+diag(ones(1,n));

