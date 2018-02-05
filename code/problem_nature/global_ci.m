function [P_CI,weights_ci] =global_ci(S1,method_,eps_)
global eps_1
eps_1 =  eps_;
nCovSamples = size(S1,3);
x0  = rand(nCovSamples,1);
x0 = x0./sum(x0);
Aeq = ones(size(x0))';
beq = 1;
lb = zeros(size(x0))';
ub = ones(size(x0))';
A =[];
b = [];
% nonlcon = [];
% f = @(x)parameterfun(x,a,b,c)
% nonlcon = @(x)circlecon(x,eps_);
nonlcon = [];


% options = optimoptions('fmincon','Display','none','Algorithm','sqp');
% options = optimoptions('fmincon','Algorithm','sqp');
% options = optimoptions('fmincon');
options = optimset('Display','none','Algorithm','sqp');
if strcmp(method_,'tr')
    x = fmincon(@cost_ci_tr,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif strcmp(method_,'det')
    x = fmincon(@cost_ci_det,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif strcmp(method_,'logdet')
    x = fmincon(@cost_ci_logdet,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end
weights_ci = x;
information_matrix_ = zeros(size(S1(:,:,1)));





for i=1:length(x)
    information_matrix_ = information_matrix_ +x(i)*inv(S1(:,:,i));
end
P_CI = inv(information_matrix_);

    function [c,ceq] = circlecon(x,eps_)
        %         global eps_1
        c = [];
        
        ceq = sum(x~=0 & x<eps_);
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
    function cost_tr = cost_ci_logdet(x)
        information_matrix = zeros(size(S1(:,:,1)));
        for i_det=1:length(x)
            information_matrix = information_matrix +x(i_det,1)*inv(S1(:,:,i_det));
        end
        cost_tr = log(det(inv(information_matrix)));
    end
    function cvx_ci_()
        cvx_begin
            variables w(2)
            % objective function is the box volume
            minimize( -log(det( w(1)*inv(S1(:,:,1)) + w(2)*inv(S1(:,:,2))  )) )
            subject to
            sum(w)==1;
            0 <= w <= 1;
        cvx_end
        disp('cvx results')
    end
end
