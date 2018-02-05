function cost_tr = cost_ci_det(x)
information_matrix = zeros(size(S1(:,:,1)));
for i_det=1:length(x)
    information_matrix = information_matrix +x(i_det,1)*inv(S1(:,:,i_det));
end
cost_tr = det(inv(information_matrix));
end