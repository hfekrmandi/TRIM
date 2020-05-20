function [A]=heat_system2(n, r)
%Boundary conditions: T(0) = 0; T(N) = 0;
N= n ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%discrete system
A=zeros(N,N);
for i=2:N-1
    A(i,i-1)=r;
    A(i,i)=1-2*r;
    A(i,i+1)=r;
end

% A(N,N) = 1 - r;
% A(N, N-1) = r;
% A(1,1)= 1-2*r;
% A(1,2) = r;

A(N,N) = 1;
A(N, N-1) = 0;
A(1,1)= 1;
A(1,2) = 0;


