function weights_ci =calc_ci_weights(S1,method_)
nCovSamples = size(S1,3);
x0  = rand(nCovSamples,1);
x0 = x0./sum(x0);
Aeq = ones(size(x0))';
beq = 1;
lb = zeros(size(x0))';
ub = ones(size(x0))';
A =[];
b = [];
nonlcon = [];
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% options = optimoptions('fmincon','Algorithm','sqp');
options = optimoptions('fmincon');
if strcmp(method_,'tr')
    x = fmincon(@cost_ci_tr,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif strcmp(method_,'det')
    x = fmincon(@cost_ci_det,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end
weights_ci = x;


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
