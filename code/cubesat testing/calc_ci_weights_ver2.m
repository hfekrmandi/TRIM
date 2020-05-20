function [weights_ci,inf_mat,inf_vect] =calc_ci_weights_ver2(S1,local_inf_vec,method_)
nCovSamples = size(S1,3);
% global opt_dist
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
if sum(weights_ci)>1
    log_message('sum of w more than 1 --> renormalizing');
    weights_ci = weights_ci./sum(weights_ci);
end

inf_mat = special_dot_sum(weights_ci,S1,1);
inf_mat = (inf_mat+inf_mat')./2;
inf_vect= special_dot_sum(weights_ci,local_inf_vec,0);
[inf_mat2 ,det_ratio]= calc_inf_ci(x);
log_message(['det_ratio = ',num2str(det_ratio)])

if max(max(inf_mat2 - inf_mat))> 0.1
    disp('OO OO')
    log_message('error in calculating information matrix ');
    
end
for i_ag=1:numel(weights_ci)
end
% if opt_dist.FLAGS.debug_CI
    %     inf_mat = special_dot_sum(weights_ci,S1,1);
    
    if ~matrix_tests(inf_mat,'is_pos')
        
        log_message('information matrix is not positive definite');
        
    end
    %     disp(['sumn of weights = ',num2str(sum(x))])
    %     log_message(['sumn of weights = ',num2str(sum(x))]);
% %     for i=1:numel(x)
% %         if ~matrix_tests( (inf_mat) - inv(S1(:,:,i)),'is_pos')
% %             log_message(['Cov mat is not bigger than local cov mat ', num2str(i)]);
% %             
% %         end
% %     end
    
% end

    function cost_tr = cost_ci_tr(x)
        
        information_matrix = zeros(size(S1(:,:,1)));
        for i_tr=1:length(x)
            information_matrix = information_matrix +x(i_tr,1)*pinv(S1(:,:,i_tr));
        end
        cost_tr = trace(inv(information_matrix));
    end
    function cost_tr = cost_ci_det(x)
        information_matrix = zeros(size(S1(:,:,1)));
        for i_det=1:length(x)
            information_matrix = information_matrix +x(i_det,1)*inv(S1(:,:,i_det));
        end
        cost_tr = -log(det(information_matrix));
        if isinf(cost_tr)
            cost_tr = log(det(inv(information_matrix)));
        end
    end

    function [information_matrix,det_ratio] = calc_inf_ci(x)
        information_matrix = zeros(size(S1(:,:,1)));
        for i_det=1:length(x)
            information_matrix = information_matrix +x(i_det,1)*inv(S1(:,:,i_det));
        end
        information_matrix = 0.5*(information_matrix+information_matrix');
        %         cost_tr = det(inv(information_matrix));
        for i_det=1:length(x)
            det_ratio(i_det) = det(information_matrix)*det(S1(:,:,i_det));
%             information_matrix = 0.5*(information_matrix + information_matrix');
%             local_cov = 0.5*((S1(:,:,i_det)) + (S1(:,:,i_det))');
%             if ~matrix_tests( (information_matrix),'is_pos')
%                 log_message('information_matrix not positive dcefinite ');
%                 
%             end
%             if ~matrix_tests( (local_cov),'is_pos')
%                 log_message([' local information_matrix not positive dcefinite ', num2str(i_det)]);
%                 
%             end
%             if ~matrix_tests( ((inv(information_matrix) - local_cov) ) ,'is_pos')
%                 log_message(['Cov mat is not bigger than local cov mat ', num2str(i_det)]);
%                 
%             end
        end
        
        
    end


end
