function C = rand_cov_mat_gen_dim(dim)

a = rand(dim,dim);
C = a'*a;

% 
% sigma_ = 1;
% for i=1:dim
%     
% sigma_ = sigma_ * (0.5 + rand(1)*[3-0.5]);
% 
% end
% C = ((1/(sigma_))*GenerateCorrelationMatrix(dim));
% end