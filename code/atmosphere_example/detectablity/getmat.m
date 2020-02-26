function [a]=getmat(alarge,n,index);
% GETMAT: Function to obtain a matrix from a matrix array
%         in which the matrices are stored column wise.
%
%         [a]=getmat(alarge,n,index);
%
%         input:
%         alarge: matrix array, a(:,i) contains the ith matrix
%         stored as 1 column i.e. [column1;column2;..;columnM].
%         n     : number of rows a.
%         index : pointer to desired matrix in matrix array.
%
%         L.G. Van Willigenburg, W.L. De Koning, 10-01-02.
%
[n1,m1]=size(alarge);m=n1/n;a=[];
if index>m1;error('Index exceeds large matrix');end
if rem(m,1)~=0;error('Large matrix has illegal size');end;
for i=1:m;a=[a alarge(1+(i-1)*n:i*n,index)];end;
