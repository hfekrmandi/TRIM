function C = rand_cov_mat_gen()
sigma_x = 0.5 + rand(1)*[3-0.5];
sigma_y = 0.5 + rand(1)*[3-0.5];
sigma_z = 0.5 + rand(1)*[3-0.5];
C = ((1/(sigma_x*sigma_y*sigma_z))*GenerateCorrelationMatrix(3));
end