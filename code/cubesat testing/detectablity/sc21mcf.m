function [s1,i0]=sc21mcf(x)
% SC21MCF : Function obtains scaling factors of each matrix column
%
%           [s1,i0]=sc21mcf(x)
%
%           x   :  input matrix
%           s1  :  1./max(abs(x)) with Inf's replaced with ones
%           i0  :  indices of totally zero columns
%
% GvW 14-02-2010

s1=max(abs(x)); i0=(s1==0); s1(i0)=1; s1=1./s1;