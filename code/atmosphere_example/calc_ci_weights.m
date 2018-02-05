function [weights_ci,inf_mat] =calc_ci_weights(S1,method_)
nCovSamples = size(S1,3);
global opt_dist
x0  = rand(nCovSamples,1);
x0 = x0./sum(x0);
Aeq = ones(size(x0))';
beq = 1;
lb = zeros(size(x0))';
ub = ones(size(x0))';
A =[];
b = [];
nonlcon = [];
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
% options = optimoptions('fmincon','Algorithm','sqp');
% options = optimoptions('fmincon');
% if strcmp(method_,'tr')
%     x = fmincon(@cost_ci_tr,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% elseif strcmp(method_,'det')
%     x = fmincon(@cost_ci_det,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% end

if strcmp(method_,'tr')
    x = fmincon(@cost_ci_tr,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif strcmp(method_,'det')
    x = fmincon(@cost_ci_det,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end


weights_ci = x;

    inf_mat = special_dot_sum(weights_ci,S1,1);

if opt_dist.FLAGS.debug_CI
%     inf_mat = special_dot_sum(weights_ci,S1,1);
    
    if ~matrix_tests(A,'is_pos')
        
        log_message('information matrix is not positive definite');
        
    end
    disp(['sumn of weights = ',num2str(sum(x))])
    log_message(['sumn of weights = ',num2str(sum(x))]);
    for i=1:numel(x)
        if ~matrix_tests( inv(inf_mat) - inv(S1(:,:,i)),'is_pos')
            log_message(['Cov mat is not bigger than local cov mat ', num2str(i)]);
            
        end
    end
    
end

    function cost_tr = cost_ci_tr(x)
        
        information_matrix = zeros(size(S1(:,:,1)));
        for i_tr=1:length(x)
            information_matrix = information_matrix +x(i_tr,1)*inv(S1(:,:,i_tr));
        end
        cost_tr = trace(inv(information_matrix));
    end
    function cost_tr = cost_ci_det(x)
        information_matrix = zeros(size(S1(:,:,1)));
        for i_det=1:length(x)
            information_matrix = information_matrix +x(i_det,1)*inv(S1(:,:,i_det));
        end
        cost_tr = det(inv(information_matrix));
    end


end
