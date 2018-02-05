clear all
close all
load('case2_step_5929--17.mat')
figure
subplot(321)
h(1)=error_ellipse(A);hold on
set(h(1),'LineWidth',3)
h(2)=error_ellipse(B);hold on
h(3)=error_ellipse(C);hold on
set(h(1),'LineWidth',2,'LineStyle','--')
set(h(2),'LineWidth',2,'LineStyle','--')
set(h(3),'LineWidth',2,'LineStyle','--')

subplot(322)
h(4)=error_ellipse(S(:,:,4));hold on
h(5)=error_ellipse(S(:,:,5));hold on
h(6)=error_ellipse(S(:,:,6));hold on
set(h(4),'LineWidth',3)
set(h(5),'LineWidth',3)
set(h(6),'LineWidth',3)

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
axis equal

subplot(321)
h(10)=error_ellipse(S(:,:,10));hold on
set(h(10),'LineWidth',3,'Color',[0 0 0])

legend(h([1 2 3 10]),'A','B','C','ABC');
axis equal
subplot(322)
h(11)=error_ellipse(S(:,:,10));hold on
set(h(11),'LineWidth',3,'Color',[0 0 0])
legend(h([4 5 6 11]),'AB','BC','AC','ABC');
axis equal



subplot(324)
h(12)=error_ellipse(S(:,:,6));hold on
h(13)=error_ellipse(S(:,:,2));hold on
h(14)=error_ellipse(S(:,:,10));hold on
h(15)=error_ellipse(S(:,:,8));hold on
set(h(12),'LineWidth',3)
set(h(13),'LineWidth',3)
set(h(15),'LineWidth',3)
set(h(14),'LineWidth',3,'Color',[0 0 0])
axis equal
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
axis equal
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
axis equal
legend(h([12 13 15 14]),'BC','A','BC-A','ABC');


% 
% set(h(1),'LineWidth',2,'LineStyle','--')
% set(h(2),'LineWidth',2,'LineStyle','--')
% set(h(3),'LineWidth',2,'LineStyle','--')
% set(h(4),'LineWidth',3)
% set(h(5),'LineWidth',3)
% set(h(6),'LineWidth',3)
% set(h(7),'LineWidth',3)
% set(h(8),'LineWidth',3)
% set(h(9),'LineWidth',3,'Color',[0 0 0])



% h(8)=error_ellipse(D12);hold on
% h(9)=error_ellipse(D13);hold on
% 
% h(10)=error_ellipse(D14);hold on
% h(11)=error_ellipse(D15);hold on
% h(12)=error_ellipse(Ddif);hold on

% 
% 
% legend(h,'A','B','C','D1 = CI(A,CI(B,C)); ',...
% 'D2 = CI(CI(A,B),C);', ...
% 'D3 = CI(CI(A,C),B);',...
% 'D4 = CI(B,CI(A,C));',...
% 'D5 = CI(CI(B,C),A);',...
% 'P_IC')

% t.String = ['D1 = CI(A,CI(B,C)); ',...
% 'D2 = CI(CI(A,B),C);', ...
% 'D3 = CI(CI(A,C),B);',...
% 'D4 = CI(B,CI(A,C));',...
% 'D5 = CI(CI(B,C),A);',...
% 'D12 = (D1-D2);',...
% 'D13 = (D1-D3);',...
% 'D14 = (D1-D4);',...
% 'D15 = (D1-D5);',...
% 'Ddif = D1 - P_{CI} ;'];