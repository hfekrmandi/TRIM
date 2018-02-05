function analysis_strict_performance()
I1 = sym('I1',[2 2],'real');
I1=triu(I1,1)+triu(I1,0).' ;
Y1 = sym('Y1',[2 2],'real');
Y1=triu(Y1,1)+triu(Y1,0).' ;


I2 = sym('I2',[2 2],'real');
I2=triu(I2,1)+triu(I2,0).' ;
Y2 = sym('Y2',[2 2],'real');
Y2=triu(Y2,1)+triu(Y2,0).' ;
[w1,C1] = symb_ci(Y1,Y2)
[w2,C2] = symb_ci(Y1+I1,Y2+I2)


diff_ = w1*Y1+(1-w1)*Y2 -(w21*Y1+(1-w2)*Y2)

end
function [w,C] = symb_ci(A,B)
[V_a,D_a] = eig(A);
T1 = [1/sqrt(D_a(1,1)) 0 ; 0 1/sqrt(D_a(2,2))]*V_a.';
[V_b,D_b] = eig(T1*B*T1.');
d1_hat = 1/(D_b(1,1));
d2_hat = 1/(D_b(2,2));

d1_til = d1_hat/(1-d1_hat);
d2_til = d2_hat/(1-d2_hat);
w = (-1/2)*(d1_til + d2_til);
% if w<0 || w>1
%     if det(A)<=det(B)
%         C = A;
%     else
%         C = B;
%     end
% else
%     C = w*A + (1-w)*B;
% end
    C = w*A + (1-w)*B;


end