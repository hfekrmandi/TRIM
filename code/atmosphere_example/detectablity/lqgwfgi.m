function [WfgI,meth,mt]=lqgwfgi(str)
% LQGWFGI: Function that determines the LQG design weighing factors and the method
%          meth=0: redundancy incorporated
%          mt    : 1,2,3 method of approximating indefinite quadratic forms see cin2nn
%
% Programmer GvW/WdK 21-5-2004
global W_f W_g W_I me_th m_t;

if nargin~=1; error(' One input argument required'); end;
if ~isstr(str); error(' Input should be a string'); end;

if str=='f'; WfgI=W_f;
elseif str=='g'; WfgI=W_g;
elseif str=='phi'; WfgI=W_I;
else
  error(' Input string should be either f,g, or phi')
end
meth=me_th; mt=m_t;