function [error_norm_1,error_norm_2] = counter_exampole_association
Aeq = [];
beq = [];

%% some optimal values
% % % x =
% % %     1.0000    0.1000    0.5000    0.0500    0.1000    1.0000
% % %     0.9683    0.1563    0.5717    0.0738    0.1238    0.7028
% % %     0.9802    0.2087    0.6908    0.1029    0.1454    0.7080
% % %     0.9497    0.1804    0.6987    0.1019    0.1646    0.7282
% % %     1.4138    0.3526    1.9937    0.1364    0.3487    2.0000
% % %     0.6230    0.0961    1.8130    0.0100    0.0100    2.0000


%%


% x0 = [1 0.1 0.5 0.05 0.1 1];
% x0 = [1.4138    0.3526    1.9937    0.1364    0.3487    2.0000];
% lb = 0.01*ones(size(x0))';
% ub = 2*ones(size(x0))';
% A =[];
% b = [];
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
% for i=1:10
% x = fmincon(@counter_example_cost1,x0,A,b,Aeq,beq,lb,ub,@nonlcon_1,options);
% x0 =x;
% end
x = [0.6230    0.0961    1.8130    0.0100    0.0100    2.0000];
% cost = counter_example_cost1(x)
close all
clc
%%
A = [x(1) 0;0 x(2)];
B = [x(3) 0;0 x(4)];
C = [x(5) 0;0 x(6)];


S(:,:,1) = A;
S(:,:,2) = B;
S(:,:,3) = C;
S(:,:,4) = global_ci(S(:,:,[1 2]),'det',eps); % A B
S(:,:,5) = global_ci(S(:,:,[2 3]),'det',eps); % B C
S(:,:,6) = global_ci(S(:,:,[1 3]),'det',eps); % A C

S(:,:,7) = global_ci(S(:,:,[3 4]),'det',eps); % A B - C
S(:,:,8) = global_ci(S(:,:,[2 6]),'det',eps); % A C - B
S(:,:,9) = global_ci(S(:,:,[1 5]),'det',eps); % B C - A

S(:,:,10) = global_ci(S(:,:,[1 2 3 ]),'det',eps); % A B C


S78 = S(:,:,7)-S(:,:,8);
S89 = S(:,:,8)-S(:,:,9);
S79 = S(:,:,7)-S(:,:,9);
S107 = S(:,:,10)-S(:,:,7);

error_norm_1 = norm(S78) + norm(S89) + norm(S79) + norm(S107);

figure
subplot(211)
h(1)=error_ellipse(A);hold on
set(h(1),'LineWidth',3)
h(2)=error_ellipse(B);hold on
h(3)=error_ellipse(C);hold on
set(h(1),'LineWidth',1.5,'LineStyle','--','Color',[0 0 0])
set(h(2),'LineWidth',1.5,'LineStyle',':','Color',[0 0 0])
set(h(3),'LineWidth',1.5,'LineStyle','-.','Color',[0 0 0])
h(10)=error_ellipse(S(:,:,10));hold on
set(h(10),'LineWidth',1.5,'Color',[0 0 0])
h_legend = legend(h([1 2 3 10]),'A','B','C','$\mathbf{CI}(A,B,C)$');
set(h_legend,'Interpreter','latex')

xlim([-1.7 1.7]); ylim([-1.5 1.5]); axis equal ;
set(gca,'fontsize', 14)
subplot(212)
h(7)=error_ellipse(S(:,:,7));hold on
h(8)=error_ellipse(S(:,:,8));hold on
h(9)=error_ellipse(S(:,:,9));hold on
h(12)=error_ellipse(S(:,:,10));hold on
set(h(7),'LineWidth',1.5,'LineStyle','--','Color',[0 0 0])
set(h(8),'LineWidth',1.5,'LineStyle',':','Color',[0 0 0])
set(h(9),'LineWidth',1.5,'LineStyle','-.','Color',[0 0 0])
set(h(12),'LineWidth',1.5,'Color',[0 0 0])
h(10)=error_ellipse(S(:,:,10));hold on
set(h(10),'LineWidth',1.5,'Color',[0 0 0])
h_legend = legend(h([7 8 9 12]),'$\mathbf{CI}(\mathbf{CI}(A,B),C)$', '$\mathbf{CI}(\mathbf{CI}(A,C),B)$', '$\mathbf{CI}(\mathbf{CI}(B,C),A)$', '$\mathbf{CI}(A,B,C)$');
set(h_legend,'Interpreter','latex')
xlim([-.5 .5]); ylim([-.5 .5]); axis equal ;
set(gca,'fontsize', 14)




%%
figure
subplot(321)
h(1)=error_ellipse(A);hold on
set(h(1),'LineWidth',3)
h(2)=error_ellipse(B);hold on
h(3)=error_ellipse(C);hold on
set(h(1),'LineWidth',2,'LineStyle','--')
set(h(2),'LineWidth',2,'LineStyle','--')
set(h(3),'LineWidth',2,'LineStyle','--')

xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis equal
subplot(322)
h(4)=error_ellipse(S(:,:,4));hold on
h(5)=error_ellipse(S(:,:,5));hold on
h(6)=error_ellipse(S(:,:,6));hold on
set(h(4),'LineWidth',3)
set(h(5),'LineWidth',3)
set(h(6),'LineWidth',3)
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;

subplot(323)
h(7)=error_ellipse(S(:,:,7));hold on
h(8)=error_ellipse(S(:,:,8));hold on
h(9)=error_ellipse(S(:,:,9));hold on
h(12)=error_ellipse(S(:,:,10));hold on
set(h(7),'LineWidth',3)
set(h(8),'LineWidth',3)
set(h(9),'LineWidth',3)
set(h(12),'LineWidth',3,'Color',[0 0 0])

legend(h([7 8 9 12]),'AB-C','AC-B','BC-A','ABC');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;

subplot(321)
h(10)=error_ellipse(S(:,:,10));hold on
set(h(10),'LineWidth',3,'Color',[0 0 0])

legend(h([1 2 3 10]),'A','B','C','ABC');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;
subplot(322)
h(11)=error_ellipse(S(:,:,10));hold on
set(h(11),'LineWidth',3,'Color',[0 0 0])
legend(h([4 5 6 11]),'AB','BC','AC','ABC');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;



subplot(324)
h(12)=error_ellipse(S(:,:,6));hold on
h(13)=error_ellipse(S(:,:,2));hold on
h(14)=error_ellipse(S(:,:,10));hold on
h(15)=error_ellipse(S(:,:,8));hold on
set(h(12),'LineWidth',3)
set(h(13),'LineWidth',3)
set(h(15),'LineWidth',3)
set(h(14),'LineWidth',3,'Color',[0 0 0])
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;
legend(h([12 13 15 14]),'AC','B','AC-B','ABC');


subplot(325)
h(12)=error_ellipse(S(:,:,4));hold on
h(13)=error_ellipse(S(:,:,3));hold on
h(14)=error_ellipse(S(:,:,10));hold on
h(15)=error_ellipse(S(:,:,7));hold on
set(h(12),'LineWidth',3)
set(h(13),'LineWidth',3)
set(h(15),'LineWidth',3)
set(h(14),'LineWidth',3,'Color',[0 0 0])
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;
legend(h([12 13 15 14]),'AB','C','AB-C','ABC');



subplot(326)
h(12)=error_ellipse(S(:,:,5));hold on
h(13)=error_ellipse(S(:,:,1));hold on
h(14)=error_ellipse(S(:,:,10));hold on
h(15)=error_ellipse(S(:,:,9));hold on
set(h(12),'LineWidth',3)
set(h(13),'LineWidth',3)
set(h(15),'LineWidth',3)
set(h(14),'LineWidth',3,'Color',[0 0 0])
xlim([-1.5 1.5]); ylim([-1.5 1.5]); axis equal ;
legend(h([12 13 15 14]),'BC','A','BC-A','ABC');


end
%%
function cost = counter_example_cost1(x)
A = [x(1) 0;0 x(2)];
B = [x(3) 0;0 x(4)];
C = [x(5) 0;0 x(6)];

S(:,:,1) = A;
S(:,:,2) = B;
S(:,:,3) = C;
S(:,:,4) = global_ci(S(:,:,[1 2]),'det',eps); % A B
S(:,:,5) = global_ci(S(:,:,[2 3]),'det',eps); % B C
S(:,:,6) = global_ci(S(:,:,[1 3]),'det',eps); % A C

S(:,:,7) = global_ci(S(:,:,[3 4]),'det',eps); % A B - C
S(:,:,8) = global_ci(S(:,:,[2 6]),'det',eps); % A C - B
S(:,:,9) = global_ci(S(:,:,[1 5]),'det',eps); % B C - A

S(:,:,10) = global_ci(S(:,:,[1 2 3 ]),'det',eps); % A B C

S108 = (sqrt(prod(diag(S(:,:,8)))) - sqrt(prod(diag(S(:,:,10)))))  / sqrt(prod(diag(S(:,:,10))));

cost = -norm(S108)

end
function [c,ceq]=nonlcon_1(x)
c = -[x(1)*x(2);...
    x(3)*x(4);...
    x(5)*x(6)];
ceq = [];
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