function x = counter_example_ver2()
dbclear if error
% dbstop in error
% % % syms a1 real
% % % syms a2 real
% % % syms b1 real
% % % syms b2 real
% % % syms c1 real
% % % syms c2 real
% % % syms d1 real
% % % syms d2 real
% % % syms w1 real
% % % syms w2 real
% % %
% % % A = [a1,0;0 a2];
% % % B = [b1,0;0 b2];
% % % C = [c1,0;0 c2];
% % % D = [d1,0;0 d2];

% x0=[1 0.1 0.2 1 0.1 0.1 0.01 0.01];
x0 = [ 1.8753    2.0000    1.3876    1.3876    2.0000    2.0000    2.0000    2.0000];

% options = optimoptions('fmincon','Algorithm','sqp');
% options = optimoptions('fmincon');
% if strcmp(method_,'tr')
%     x = fmincon(@cost_ci_tr,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% elseif strcmp(method_,'det')
%     x = fmincon(@cost_ci_det,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% end
Aeq = [];
beq = [];
lb = 0.01*ones(size(x0))';
ub = [];2*ones(size(x0))';
A =[];
b = [];
nonlcon = [];
options = optimoptions('fmincon','Display','none','Algorithm','interior-point');
x = fmincon(@counter_example_cost,x0,A,b,Aeq,beq,lb,ub,@nonlcon_,options);
% x = [ 1.8753    2.0000    1.3876    1.3876    2.0000    2.0000    2.0000    2.0000];
S1(:,:,1) = [x(1) 0 ;0 x(2)];
S1(:,:,2) = [x(3) 0 ;0 x(4)];
S1_2(:,:,1) = [x(1)+x(5) 0;0 x(2)+x(6) ];
S1_2(:,:,2) = [x(3)+x(7) 0;0 x(4)+x(8) ];
method_ = 'det'
[weights_hyb,inf_mat_hyb] =ci_weights_nature(S1,method_);
[weights_ci,inf_mat_ci] =ci_weights_nature(S1_2,method_);
delta_I1 = [x(5) 0;0 x(6)];
delta_I2 = [x(7) 0;0 x(8)];
cost_tr = det(inf_mat_hyb +delta_I1+ delta_I2) - det(inf_mat_ci)  ;


    function cost_tr = counter_example_cost(x)
        S1(:,:,1) = [x(1) 0 ;0 x(2)];
        S1(:,:,2) = [x(3) 0 ;0 x(4)];
        S1_2(:,:,1) = [x(1)+x(5) 0;0 x(2)+x(6) ];
        S1_2(:,:,2) = [x(3)+x(7) 0;0 x(4)+x(8) ];
        method_ = 'det';
        [weights_hyb,inf_mat_hyb] =ci_weights_nature(S1,method_);
        [weights_ci,inf_mat_ci] =ci_weights_nature(S1_2,method_);
        delta_I1 = [x(5) 0;0 x(6)];
        delta_I2 = [x(7) 0;0 x(8)];
        cost_tr = det(inf_mat_hyb +delta_I1+ delta_I2) - det(inf_mat_ci)  
        
        
    end
    function [c,ceq]=nonlcon_(x)
        I1 = [x(1) 0 ;0 x(2)];
        I2 = [x(3) 0 ;0 x(4)];
        I3 = [x(5) 0;0 x(6) ];
        I4 = [x(7) 0;0 x(8) ];
        
        c(1) = -det(I1)+0.01;
        c(2) = -det(I2)+0.01;
        c(3) = -det(I3)+0.01;
        c(4) = -det(I4)+0.01;
        c(5) = -det(I4-I3)+0.01;
        c(6) = -det(I4-I2)+0.01;
        c(7) = -det(I4-I1)+0.01;
        c(8) = -det(I3-I2)+0.01;
        c(9) = -det(I3-I1)+0.01;
        c(10) = -det(I2-I1)+0.01;
        c(11) = (det(inf_mat_hyb +delta_I1+ delta_I2) - det(inf_mat_ci))+0.01; 
        ceq = [];
    end

end