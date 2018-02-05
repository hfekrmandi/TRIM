
NUM_SYS             = 30;   %   DOF of system
L = 1;                      %   Length of the rod
alf = 4.2*1e-5;             %   thermal diffusivity
xlim = [0, L];
dx = (xlim(2) - xlim(1))/ (NUM_SYS); 
nNodes = 10;
x0 = xlim(1) + [0 : nNodes] * dx;  % nodes 
dt = 4;                         % time step
r = alf * dt/(dx^2);            % unstable if r > 0.5
[A]=heat_system2(NUM_SYS ,r); 
p = eye(size(A,1));
for i=1:100
    p = A'*p*A;
end