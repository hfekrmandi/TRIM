covariance_ = [0.8 0.2;0.2 0.3];
mu =[1;2];
close all
clear all
close all
color = [1 0 0];
alpha_value = 0.2
mu = [0;0];
theta = pi/4
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
cov_init = [1.5 0;0 .01];
% figure
cov_init3 = [1 0;0 .01];
S1(:,:,1) = cov_init3;
h1=error_ellipse(S1(:,:,1),mu) ; hold on
set(h1,'Linewidth',3)
color = [0 1 0];
cov_init2 = [0.1 0;0 .01];;

S1(:,:,2) = R^2*cov_init*R^2';
h2= error_ellipse(S1(:,:,2),mu)
set(h2,'Linewidth',3)
color = [0 0 1];

S1(:,:,3) = R*cov_init2*(R');
h3= error_ellipse(S1(:,:,3),mu)
set(h3,'Linewidth',3)

S1(:,:,4) = R^3*cov_init2*(R')^3;
h4= error_ellipse(S1(:,:,4),mu)
set(h4,'Linewidth',3)


P_GCI = global_ci(S1,'logdet',eps);
h4= error_ellipse(P_GCI,mu)
set(h4,'Color',0*[1 1 1])
set(h4,'Linewidth',3)
P_MH = inv(inv(S1(:,:,1)) + inv(S1(:,:,2)) +  inv(S1(:,:,3)) + inv(S1(:,:,4)) )
h5= error_ellipse(P_MH,mu)
set(h5,'Color',1*[1 0 1])
set(h5,'Linewidth',3)
set(h5,'LineStyle','-.')

axis equal

% for i_agent=1:3
% %     S1(:,:,i_agent) = rand_cov_mat_gen_dim(dim);
% patch_ellipse(S1(:,:,i_agent),mu,color,alpha_value)
% end
% P_GCI = global_ci(S1,'det',eps_);
% 
% patch_ellipse(covariance_,mu,color,alpha_value)