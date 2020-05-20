function[al,si]=sric1(phi,gam,qi,ri,am,si)
%SRIC1
%
%    [al,si]=sric1(phi,gam,qi,ri,am,si)
%
%    function realizes one recursion of the discrete riccati equation
%    resulting from a regulator problem with quadratic cost including
%    a crossterm AM (see edorti.m and edortv.m)
%
%    PHI,GAM   : equivalent discrete system matrices
%    QI,RI,AM  : equivalent matrices quadratic cost (including crossterm AM)
%    SI        : solution to riccati equation
%    AL        : optimal feedback

gsgr=gam'*si*gam+ri; gspm=gam'*si*phi+am'; psp=phi'*si*phi;
al=pinv(gsgr)*gspm; sh=qi+psp-al'*gspm; %sh=qi+psp-al'*gsgr*al;
si=0.5*(sh+sh'); % Recover symmetry
%sh=msqrt(sh); si=sh*sh'; % Recover symmetry