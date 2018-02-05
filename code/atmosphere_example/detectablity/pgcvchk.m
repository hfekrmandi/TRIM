% PGCVCHK : Check time-varying dimensions of equivalent discrete-time
%           LQ control problem in the case of aynchronous sampling
%
% GvW: 13-10-2003

  function [ni,mi,ni1,li]=pgcvchk(phi,gam,qi,ri,mc,vi,c,w,me)
  
  if ~(nargin==6 || nargin==9); error(' 6 or 9 input arguments required'); end;
  [ni1,ni]=size(phi); [nh,mi]=size(gam);
  if nh~=ni1 && nh~=0; error(' Dimensions gam inconsistent with phi'); end;
  [nh,mh]=size(qi);
  if nh~=ni || mh~=ni; error(' Dimensions qi inconsistent with phi'); end;
  [nh,mh]=size(ri);
  if nh~=mi || mh~=mi; error(' Dimensions ri inconsistent with gam'); end;
  [nh,mh]=size(mc);
  if (nh~=ni && nh~=0) || mh~=mi; error(' Dimensions mc inconsistent with phi,gam'); end;
  [nh,mh]=size(vi);
  if nh~=ni1 || mh~=ni1; error(' Dimensions vi inconsistent with phi'); end;
  
  if nargin > 6
     [li,mh]=size(c);
     if min(li,mh)~=0 && mh~=ni; error(' Dimensions c inconsistent with p'); end;
     [nh,mh]=size(w);
     if min(nh,mh)~=0 && nh~=li || mh~=li; error(' Dimensions w inconsistent with c'); end;
     [nh,mh]=size(me);
     if (nh~=ni1 && li~=0) || mh~=li ; error(' Dimensions me inconsistent with p, c'); end;
  end
