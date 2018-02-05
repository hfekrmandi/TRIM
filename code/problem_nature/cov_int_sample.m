function [error_norm_1,error_norm_2] = cov_int_sample(dim)
persistent case_no_1  case_no_2 i_iter
if isempty(case_no_2)
   case_no_2 = 1;
end

if isempty(i_iter)
   i_iter = 0;
end

if isempty(case_no_1)
   case_no_1 = 1;
end

i_iter = i_iter+1;


A = rand_cov_mat_gen_dim(dim);
B = rand_cov_mat_gen_dim(dim);
C = rand_cov_mat_gen_dim(dim);
S(:,:,1) = A;
S(:,:,2) = B;
S(:,:,3) = C;
S(:,:,4) = global_ci(S(:,:,[1 2]),'det'); % A B
S(:,:,5) = global_ci(S(:,:,[2 3]),'det'); % B C 
S(:,:,6) = global_ci(S(:,:,[1 3]),'det'); % A C

S(:,:,7) = global_ci(S(:,:,[3 4]),'det'); % A B - C
S(:,:,8) = global_ci(S(:,:,[2 6]),'det'); % A C - B
S(:,:,9) = global_ci(S(:,:,[1 5]),'det'); % B C - A

S(:,:,10) = global_ci(S(:,:,[1 2 3 ]),'det'); % A B C


S78 = S(:,:,7)-S(:,:,8);
S89 = S(:,:,8)-S(:,:,9);
S79 = S(:,:,7)-S(:,:,9);
S107 = S(:,:,10)-S(:,:,7);

error_norm_1 = norm(S78) + norm(S89) + norm(S79) + norm(S107); 
    





% 
% 
% D1 = cov_int(A,cov_int(B,C));
% D2 = cov_int(cov_int(A,B),C);
% D3 = cov_int(cov_int(A,C),B);
% D4 = cov_int(B,cov_int(A,C));
% D5 = cov_int(cov_int(B,C),A);
% 
% 
% D12 = (D1-D2);
% D13 = (D1-D3);
% D14 = (D1-D4);
% D15 = (D1-D5);
% D_dif = D1 - S(:,:,10) ;
% 
% error_norm_2 = norm(D12) + norm(D13) + norm(D14) + norm(D15) + norm(D_dif);
error_norm_2 =error_norm_1;
if error_norm_2>0.2
    disp(error_norm_2)
    save(['case1_step_',num2str(i_iter),'--',num2str(case_no_1)])
    case_no_1 = case_no_1 + 1;
end

if error_norm_1>0.2
    disp(error_norm_1)
    save(['case2_step_',num2str(i_iter),'--',num2str(case_no_2)])
    case_no_2 = case_no_2 + 1;

end


end
function C = cov_int(A,B)
[V_a,D_a] = eig(A);
T1 = [1/sqrt(D_a(1,1)) 0 ; 0 1/sqrt(D_a(2,2))]*V_a.';
[V_b,D_b] = eig(T1*B*T1.');
d1_hat = 1/(D_b(1,1));
d2_hat = 1/(D_b(2,2));

d1_til = d1_hat/(1-d1_hat);
d2_til = d2_hat/(1-d2_hat);
w = (-1/2)*(d1_til + d2_til);
if w<0 || w>1
    if det(A)<=det(B)
        C = A;
    else
        C = B;
    end
else
    C = inv(w*inv(A) + (1-w)*inv((B)));
end

end