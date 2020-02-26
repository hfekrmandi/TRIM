function [Ar,Br,Cr] = create_reduced_atmospheric_model(dim_sys)
% This function takes a dim_sys which is supposed to be less than 810
% (I recommend below 240) and first builds an atmospheric system and 
% then reduces it to a system with A of dim_sys * dim_sys


%Source and receptors are not on the boundaries (except z = 0)
%10^5 system,
%1. Full order system: A(SR, SC, S), B
%2. B: 10 point sources
%3. C: You can define the measurement yourself, C=  zeros(NUM_OUT,
%NUM_SYS); place the sensors you want
%4. ROM : order of the reduced model can be changed from 100 - 380, through
%NUM_NON, NUM_NON = 150 is enough for most cases, if you want higher
%accuracy, try NUM_NON  = 380.
%5. There are two set of ROM you can use,
%either you (Ar, Br, Cr, UF, VF), which you might have complexed valued
%system, or you can you (At, Bt, Ct, Tlt, Trt).
%5. Actual field is generated through X at the end of the code, change f if
%you want different constant forcing, the process noise is eta, and is
%added in all the field.
%6. Don't forget to use the transformation when you do the KF on the ROM
%%
%SET PARAMETERS : Parameters can be changed
%Computational domain
xlim                            =   [ 0, 2000];  %[600m - 3000m]
ylim                            =   [-100, 400];
zlim                            =   [0, 50];
%Wind Profile parameters
Uwind                           =   4;    % wind speed (m/s) between 1-5 (m/s)
alp                             =   0/180 * pi;
beta                            =   0/180 * pi;
%Diffusion parameters
ay                              =   0.08 *0.1;
by                              =   0.0001 * 0.1;
az                              =   0.06 * 0.1;
bz                              =   0.0015 * 0.1;
%discretization dimension
nx                              =   100;
ny                              =   100;
nz                              =   10;
dt                              =   1;

%%
%U = (U cos(alp)cos(beta), U cos(alpha) sin(beta), U sin(alp));
Ux = Uwind * cos(alp)*cos(beta);
Uy = Uwind * cos(alp)*sin(beta);
Uz = Uwind * sin(alp);
%Finite difference Scheme parameters
dx =( xlim(2) - xlim(1))/nx;
dy = (ylim(2) - ylim(1))/ny;
dz = (zlim(2) - zlim(1))/nz;
x0 = xlim(1) + [0 : nx]*dx;           % distance along wind direction (m)
y0 = ylim(1) + [0 : ny]*dy;           % cross-wind distance (m)
z0 = zlim(1) + [0 : nz]*dz;
[xmesh, ymesh] = meshgrid(x0,y0);
Ky = zeros(nx+1,1);
Kz = zeros(nx+1, 1);
for i = 1 : nx+1
    Ky(i,1) = Uwind * ay^2/2 * (2 * x0(i) + by *(x0(i))^2)/(1+by*x0(i))^2;
    Kz(i,1) = Uwind * az^2/2 * (2 * x0(i) + bz *(x0(i))^2)/(1+bz*x0(i))^2;
end

sy = Ky * dt/(dy^2);
sz = Kz * dt/(dz^2);
cx = Ux * dt/dx;
cy = Uy * dt/(dy);
cz = Uz * dt/(dz);
%Stable

NUM_SYS = (ny-1)*(nx-1)*(nz);
%Primal system
[SR, S, SC] = atmo3dlax2 (Uwind, alp, beta, nx, ny, nz, xlim, ylim, zlim, Ky, Kz, dt);
%Adjoint system : Sad = S;
SRad = SC;
SCad = SR;
%%
%Source and receptor locations and discretization
% Stack emission source data:
source.n = 10;                         % # of sources
source.x = [280, 300, 900, 1100, 500,1300,500,1700,1400,700];     % x-location (m)
source.y = [ 75, 205, 25,  185,  250,230, 330,170, 240, 200];     % y-location (m)

source.z = [ 15,  35,  15,   15, 10, 10, 10, 10, 15, 20];     % height (m)
%source.label=[' S1'; ' S2'; ' S3'; ' S4'; 'S5';'S6';'S7';'S8';'S9';'S10'];


tpy2kgps = 1.0 / 31536;               % conversion factor (tonne/yr to kg/s)
source.Q = [35, 80, 5, 5, 10, 5, 5, 10, 5, 10] * tpy2kgps ; % emission rate (kg/s)
%source.Q = 1;
% Stack discrete source data: Transfer from continuous to discrete
dissource.n = source.n;
dissource.x = floor((source.x - xlim(1))/dx)+1;
dissource.y = floor((source.y - ylim(1))/dy)+1;
dissource.z = floor((source.z - zlim(1))/dz)+1;
%%
B = zeros(NUM_SYS, source.n);
for ns = 1 : source.n
    B((nx-1)*(ny-1)*(dissource.z(1,ns)-1)+(nx-1)*(dissource.y(1,ns)-2)+dissource.x(1,ns)-1, ns) = 1;
end


NUM_IN = source.n;
NUM_OUT = NUM_SYS;
%C = eye(NUM_SYS);
%D = zeros(NUM_OUT, NUM_IN);
%sys = ss(A, B, C, D, 1);
Ft = (1 - (source.x - x0(dissource.x))/dx).*(1 -(source.y - y0(dissource.y))/dy) .* (1- (source.z - z0(dissource.z))/dz) .* source.Q *dt/(dx*dy*dz);
C = zeros(dissource.n, NUM_SYS);
for ns = 1 : dissource.n
    C(ns,(nx-1)*(ny-1) * (dissource.z(1,ns)-1) + (nx-1) * (dissource.y(1,ns) -2) + dissource.x(1,ns)-1) = 1;
end

%%
%Comparision with the puff plume solution
%%
%RPOD

%NUM_IN                  = NUM_SYS;
%NUM_OUT                 = recept.n;

%start with the full field measurements
%Collect snapshots at seperated timestep.
%F = (randn(source.n, TIME_STEP * NUM_SNAP));
%%


TIME_STEP = 10;
NUM_ROUND = 401;

tic
%Xt = zeros(NUM_SYS, NUM_ROUND * TIME_STEP);
%for i = 1 : NUM_ROUND*TIME_STEP
%    Xt(:,i+1) = A*Xt(:,i) + B*randn(NUM_IN,1);
%end
%Xs = Xt(:, TIME_STEP:TIME_STEP:300 *TIME_STEP);
Xt = zeros(NUM_SYS, 1);
for iall = 1 : NUM_ROUND
    X = zeros(NUM_SYS,TIME_STEP);
    X(:,1) = Xt;
    for k= 1 : TIME_STEP - 1
        for i = 1  : length(SR)
            X(SR(i),k+1) = X(SR(i),k+1)+S(i)*X(SC(i),k);
        end
        X(:,k+1) = X(:,k+1) + B*randn(NUM_IN,1);
    end
    Xt = X(:,end);
    Xs(:, iall) = X(:,end); %Pick (X(TIME_STEP, 2*TIME_STEP,..., NUM_ROUND*TIME_STEP))
end
rpod_genx = toc
tic
%TIME_STEP = 20;
%NUM_SNAP = 100;
Yt = zeros(NUM_SYS, 1);
%Yp = zeros(N, NUM_ROUND);
%F = randn(disrecept.n, TIME_STEP * NUM_SNAP);
%F = abs(randn(N, TIME_STEP*NUM_SNAP));
for iall = 1 : NUM_ROUND
    Y = zeros(NUM_SYS,TIME_STEP);
    Y(:,1) = Yt;
    for k= 1 : TIME_STEP - 1
        for i = 1  : length(SC)
            Y(SC(i),k+1) = Y(SC(i),k+1)+S(i)*Y(SR(i),k);
            
        end
        %Y(:,k+1) = Y(:,k+1) + F(:, TIME_STEP*(iall-1)+k);
        Y(:,k+1) = Y(:,k+1) + (randn(NUM_SYS,1));
    end
    Yt = Y(:,end);
    Ys(:, iall) = Y(:,end);
    %Yp(:, iall) = Y(:,end);
end
%Yt = zeros(NUM_SYS, NUM_ROUND * TIME_STEP);
%for i = 1 : NUM_ROUND*TIME_STEP
%    Yt(:,i+1) = A' * Yt(:,i) + randn(NUM_SYS,1);
%end
%Ys = Yt(:, TIME_STEP: TIME_STEP: 300 *TIME_STEP);
rpod_geny = toc
% save x1data Xs Ys
clear Xt Yt X Y
Xp = Xs(:, 1:end);
Yp = Ys(:, 1:end);

tic
Ht =Yp' * Xp;
rpod_hankel = toc
%%
tic
[Ut, St, Vt] = svd(Ht);
rpod_svd = toc
NUM_NON = 80;
Unt = Ut(:, 1:NUM_NON);
Snt = St(1:NUM_NON, 1:NUM_NON);
Vnt = Vt(:,1:NUM_NON);

Trt= Xp* Vnt* Snt^(-1/2);
Tlt = Snt^(-1/2) * Unt' * Yp';

%clear X  Y  Xt Yt

Temp = zeros(NUM_SYS,NUM_NON);
for j = 1 : NUM_NON
    for i = 1 : length(SR)
        Temp (SR(i), j) = Temp (SR(i),j) + S(i) * Trt(SC(i),j);
    end
end

At=  Tlt*Temp;
Bt = Tlt * B;
%Ct = C * Trt;
Ct =  Trt;

% clear Temp
[ce,dat] = eigs(At, NUM_NON);

d3 = diag(dat);
%[ce,de] = eig(At);
VF = Trt * ce;
UF = inv(ce)* Tlt;
Ar = inv(ce)*At * ce;
Br = UF * B;
% Cr = VF;
idx_c = [24934       24935       26102       26458       28042       31824       32783       35002       36106];
Cr  = VF(idx_c',:);

% save RPODmodel Ar Br Cr UF VF

sig_w = 0.1;
TIME = 1000;
eta = sig_w * randn(NUM_SYS, TIME);
X = zeros(NUM_SYS, TIME);

for i = 1 : NUM_IN
    f(i,1) =10;
end
for k= 1 : TIME-1
    for i = 1  : length(SR)
        X(SR(i),k+1) = X(SR(i),k+1)+S(i)*X(SC(i),k);
    end
    X(:,k+1) = X(:,k+1) + B * f  + eta(:,k);
end
% save dataRPOD X sig_w f
