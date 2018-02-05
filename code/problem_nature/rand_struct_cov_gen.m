function C = rand_struct_cov_gen(dim,S,eig_)

a = rand(dim,dim).*S;
C = a'+a ;
[U,s,v] = svd(C);
C = U*diag(eig_)*U';

% 
% sigma_ = 1;
% for i=1:dim
%     
% sigma_ = sigma_ * (0.5 + rand(1)*[3-0.5]);
% 
% end
% C = ((1/(sigma_))*GenerateCorrelationMatrix(dim));
% end