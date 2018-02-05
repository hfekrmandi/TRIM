% ADJ2RC : Function that adjusts the state to include the running costs
%
%          [xout]=adj2rc(xin,dimxd)
%
%          Input:
%                xin   :  state [x;d]
%                dimxd :  dimension state without running costs
%
%          Output:
%                xout  :  [x;0;d]
%
% GvW 11-5-2005
%
function [xout]=adj2rc(xin,dimxd)

if nargin~=2; error(' 2 input arguments required'); end
if nargout~=1; error(' 1 output argument required'); end

[n,m]=size(xin);
if n~=dimxd; error(' dimension x unequal to dimxd'); end

xout=[xin;zeros(1,m)];