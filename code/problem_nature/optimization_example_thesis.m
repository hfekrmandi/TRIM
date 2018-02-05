function [P_cc,P_ab] = optimization_example_thesis(P_aa,P_bb)
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = [0,0,0];
x = fmincon(@cost_fun,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

        P_ab = [x(1) x(2);...
            x(2)  x(3)];
        P = [P_aa, P_ab;...
            P_ab',P_bb];
        I_cc = [eye(2), eye(2)]* inv(P) * [eye(2);eye(2)];
        P_cc = inv(I_cc);



    function [c,ceq] = nonlcon(x)
        P_ab = [x(1) x(2);...
            x(2)  x(3)];
        P = [P_aa, P_ab;...
            P_ab',P_bb];
 
        ceq = all(eig(P)>0) - 1;
        c = [];
    end



    function cost = cost_fun(x)
        P_ab = [x(1) x(2);...
            x(2)  x(3)];
        P = [P_aa, P_ab;...
            P_ab',P_bb];
        I_cc = [eye(2), eye(2)]* inv(P) * [eye(2);eye(2)];
        P_cc = inv(I_cc);
        cost = prod(eig(P_cc));
    end
end