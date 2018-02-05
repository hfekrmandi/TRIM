covariance_ = [0.8 0.2;0.2 0.3];
close all
clear all
close all


%% Defining covariances
mu =[0;0];
P_aa = diag([1 0.1]);
P_bb = diag([0.1 1]);

%% solving an ioptimization problem to find the extreme case
% [P_cc,P_ab] =  optimization_example_thesis(P_aa,P_bb);



%% calculate the fusion
S1(:,:,1) = P_aa;
S1(:,:,2) = P_bb;
% P_ab= [-0.3162 0 ; 0 -0.3162]
P_ab = [0.3162 0 ; 0 0.3162];

P = [P_aa, P_ab;...
     P_ab',P_bb];
I_cc = [eye(2), eye(2)]* inv(P) * [eye(2);eye(2)];


digits(4); 
ss=sym(P,'d'); 
latex(ss);
ss = sym(P);
v = vpa(ss,5); 
latex(ss)

P_GCI = global_ci(S1,'logdet',eps);
P_MH = inv(inv(S1(:,:,1)) + inv(S1(:,:,2)));
P_cc = inv(I_cc);


%% plot 1: unknown correlation solved with CI
figure
h1=error_ellipse(P_aa,mu) ; hold on
set(h1,'Linewidth',1.5)
set(h1,'LineStyle','--')
set(h1,'Color',0*[1 1 1])

h2= error_ellipse(P_bb,mu);
set(h2,'Linewidth',1.5)
set(h2,'LineStyle','--')
set(h2,'Color',0*[1 1 1])


h3= error_ellipse(P_GCI,mu);
set(h3,'Color',0*[1 1 1])
set(h3,'Linewidth',1.5)

h_legend = legend([h1 h2 h3],{'P_{aa}','P_{bb}','P_{cc} with unknown correlation P_{ab}'})
set(h_legend,'FontSize',14);
set(gca,'FontSize',16)


axis equal
set(gca,'DataAspectRatio',[1 1 1],'FontSize',16,'PlotBoxAspectRatio',...
    [434 342.3 228.2]);
% fig2texPS()
%% plot 2: no correlation 
figure
h1=error_ellipse(P_aa,mu) ; hold on
set(h1,'Linewidth',1.5)
set(h1,'LineStyle','--')
set(h1,'Color',0*[1 1 1])

h2= error_ellipse(P_bb,mu);
set(h2,'Linewidth',1.5)
set(h2,'LineStyle','--')
set(h2,'Color',0*[1 1 1])


h3= error_ellipse(P_MH,mu);
set(h3,'Color',0*[1 1 1])
set(h3,'Linewidth',1.5)

legend([h1 h2 h3],{'P_{aa}','P_{bb}','P_{cc} with zero correlation P_{ab}=0'})

axis equal
set(gca,'DataAspectRatio',[1 1 1],'FontSize',16,'PlotBoxAspectRatio',...
    [434 342.3 228.2]);
% fig2texPS()
%% plot 3: unknown correlation solved with CI
figure
h1=error_ellipse(P_aa,mu) ; hold on
set(h1,'Linewidth',1.5)
set(h1,'LineStyle','--')
set(h1,'Color',0*[1 1 1])

h2= error_ellipse(P_bb,mu);
set(h2,'Linewidth',1.5)
set(h2,'LineStyle','--')
set(h2,'Color',0*[1 1 1])

h3= error_ellipse(P_cc,mu);
set(h3,'Color',0*[1 1 1])
set(h3,'Linewidth',1.5)

legend([h1 h2 h3],{'P_{aa}','P_{bb}','P_{cc} with known correlation P_{ab}'})

axis equal
set(gca,'DataAspectRatio',[1 1 1],'FontSize',16,'PlotBoxAspectRatio',...
    [434 342.3 228.2]);
% fig2texPS()




% h5= error_ellipse(P_MH,mu)
% set(h5,'Color',1*[1 0 0])
% set(h5,'Linewidth',1.5)
% set(h5,'LineStyle','-.')
% 
% h6= error_ellipse(P_cc,mu)
% set(h6,'Color',[0 0 1])
% set(h6,'Linewidth',2.5)
% % set(h6,'LineStyle','--')
% 
% legend([h1,h6],{'a','b'})
% 

axis equal





% alpha_value = 0.2
% mu = [0;0];
% theta = pi/4

% A = rand_cov_mat_gen_dim(4);
% S = [1 0 1 1;...
%      0 1 1 1;...
%      1 1 1 0;...
%      1 1 0 1];
%  A = rand_struct_cov_gen(4,S,[1 0.1 0.1 1]);
% [Q,R] = qr(A);
% A = Q*Q';
% [U,S,V] = svd(A)
% A =V*S*V'
% cov_init = [1.5 0;0 .01];
% cov_init3 = [1 0;0 .1];
% S1(:,:,1) = A(1:2,1:2);

% h1=error_ellipse(S1(:,:,1),mu) ; hold on
% set(h1,'Linewidth',3)
% cov_init2 = [0.1 0;0 .01];
% 
% 
% S1(:,:,2) = A(3:4,3:4);
% h2= error_ellipse(S1(:,:,2),mu)
% set(h2,'Linewidth',3)
% color = [0 0 1];
% 
% 
% P_aa = S1(:,:,1);
% P_bb = S1(:,:,2);
% R = GenerateCorrelationMatrix(2) ;
% P_ab =  R;%0.09* P_aa ; 

% P_ab =  A(1:2,3:4) ; 
% 
% P = [P_aa, P_ab;...
%      P_ab',P_bb];
% I_cc1 = [eye(2), eye(2)]* inv(P) * [eye(2);eye(2)];
% I_cc2 = inv(P_aa) + (inv(P_aa)*P_ab - eye(2)) * inv(P_bb  - P_ab'* inv(P_aa) * P_ab) * (P_ab'*inv(P_aa) - eye(2));
% % P_cc2 = (P_aa) - (P_ab) * inv(P_bb) * (P_ab');
% P_cc2 = inv(I_cc2);
% 
% P_cc1 = inv(I_cc1);





