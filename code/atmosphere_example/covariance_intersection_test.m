function [est_inf1,est_inf2] = covariance_intersection_test()
adj_mat = [1 1 0;
    1 1 1;
    0 1 1] ;

init_inf = [1;1;1];

 est_inf1(:,1) = init_inf;
 est_inf2(:,1) = init_inf;



for i_iter =2:100
    est_inf1(:,i_iter) = ci_calc(est_inf1(:,i_iter - 1),adj_mat,1);
    est_inf2(:,i_iter) = ci_calc(est_inf2(:,i_iter - 1),adj_mat,2);
    
    e_diff(i_iter) = norm(est_inf1(:,i_iter)  - est_inf2(:,i_iter));
end
e_diff



end
function est_inf = ci_calc(inf_vector,adj_mat,flag_)

row_vec = 1./inf_vector';
raw_mat = repmat(row_vec,numel(inf_vector),1);
if flag_==1
w1_mat = 1./sum(row_vec).*ones(size(adj_mat));
est_inf = w1_mat*inf_vector;
end
if flag_==2
w2_mat = repmat(1./sum( raw_mat.*adj_mat,2),1,size(adj_mat,2));
est_inf = (w2_mat.*raw_mat)*inf_vector;
end
end