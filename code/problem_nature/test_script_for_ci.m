aa=  inv(S1(:,:,i)) +  inv(S1(:,:,i))';
bb = eig(2*inf_mat -ans)
figure
imagesc(1/(opt_dist.sim.obs.r_CI_var{i_agent})*(H_CI'*H_CI)+opt_dist.result.pred.Y_bar_CI(:,:,i_agent))

ss = 1/(opt_dist.sim.obs.r_CI_var{i_agent})*(H_CI'*H_CI)
max(max(ss))

information_matrix = 0.5*(information_matrix + information_matrix')
inf_ = 0.5*(inv(S1(:,:,i_det)) + inv(S1(:,:,i_det))')
matrix_tests( (information_matrix) - inf_,'is_pos')


figure
subplot(211)
imagesc(information_matrix);
subplot(212)
imagesc(local_inf);


figure
subplot(211)
imagesc(inv(information_matrix));
subplot(212)
imagesc(inv(local_inf));

% figure
% imagesc(local_cov)

%  lets try global cov intersection
for i_agent=1:9
    p_local_CI_test(:,:,i_agent) = inv(opt_dist.result.consenus{i_agent}.delta_I_CI{1});
    i_local_test(:,:,i_agent) = opt_dist.result.consenus{i_agent}.delta_I_CI{1};
end

[weights_ci_test,inf_mat_test,inf_vect_test] =calc_ci_weights_ver2(p_local_CI_test,i_local_test,'det');




% calc centralized one
delta_i_cen = zeros(size(opt_dist.result.pred.x_cen));
delta_I_cen  = zeros(size(opt_dist.result.pred.P_cen));

for i_agent = 1 : opt_dist.nAgents
    delta_i_cen = delta_i_cen  +  opt_dist.result.consenus{i_agent}.delta_i{1};
    delta_I_cen = delta_I_cen  +  opt_dist.result.consenus{i_agent}.delta_I{1};
end

% calc it different way correct(
[av_delta_I,av_delta_i] = calculate_aver_info(opt_dist.result);

updat_cen_Y = opt_dist.result.pred.Y_cen +delta_I_cen;
updat_cen_y = opt_dist.result.pred.y_cen +delta_i_cen;

close all
figure
% subplot(411)
for i=1:size(p_local_CI_test,3)
h(i) = error_ellipse(p_local_CI_test(:,:,i));hold on
set(h(i),'LineWidth',1)
% pause(2)
end

h(i+1)= error_ellipse(inv(inf_mat_test));hold on
set(h(i+1),'LineWidth',2)

h(i+2)= error_ellipse(inv(updat_cen_Y));hold on
set(h(i+2),'LineWidth',3)

h(i+3)= error_ellipse(inv(inf_mat_ci));hold on
set(h(i+3),'LineWidth',3)
% title('global CI vs centralized')
mu1 =[1;1];
mu2 =[1;1];
dist1 = BC_distance(mu1,inf_mat_test,mu2,updat_cen_Y);
dist2 = BC_distance(mu1,inf_mat_ci,mu2,updat_cen_Y);


legend([h(i+1) h(i+2) h(i+3)],'global CI','centralized','local CI one iteration')
title({['det(CI)/det(cen) = ',num2str(det(inf_mat_test)/det(updat_cen_Y)), '  BC_sit =',num2str(dist1) ],...
        ['det( local CI)/det(cen) = ',num2str(det(inf_mat_ci)/det(updat_cen_Y)),'  BC_sit =',num2str(dist2)  ]})




