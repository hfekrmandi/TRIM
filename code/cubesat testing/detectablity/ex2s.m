% EX2 : Example 2. 
%       Script which calls opconstr.m to solve the digital optimal
%       control problem specified in dopcoex.m using automatic
%       differentiation and next computes an optimal LQG compensator

load ex2
hf=figure; plot(tt,xt(:,1:end-1)); hold on; stair(tu(:,1)',tu(1:end-1,2:end)','m');
hold off; xlabel('Time'); ylabel('States');
str=['Optimal control & state trajectory min(J)=' num2str(J)];
str=[str ' tf=' num2str(tu(end,1))]; title(str);
print(hf,'Fig4.emf','-dmeta'); % Save figure
%% Compute equivalent discrete-time LQG problem data
clear; lqgcomp('tlsdyn1','tlsout1','tlsphi1','tlslqgc','dop2','montcar',3);
% Temporal stabilizability and detectability analysis
% Compute equivalent discrete-time LQG problem data
% clear; edodata('tlsdyn','tlsout','tlsphi','tlslqgc','dop1')
load edorp; % Load equivalent discrete-time LQG problem data
js=3; ks=js; tolWM=1e-6; epsrw=1e-6;
[npt,nst,ptt,stt,NW,NM,rW,rM,nAD,hfs]=tmpstdt(lqgarr,js,ks,tolWM,epsrw,epsrw);
print(hfs(1),'Fig5.emf','-dmeta');
nsd=(nst(1:end-1)-nst(2:end))./nst(2:end);
hf=figure; plot(NW(1:end-1),nsd,'o');
title('Relative one-step stabilizability measure ||S_{i}^{\epsilon}||');
xlabel('Discrete-time'); ylabel('Relative one-step stabilizability measure ||S_{i}^{\epsilon}||')
print(hf,'Fig6.emf','-dmeta');
% Display table
disp('     i rank(W) rank(M) nAD'); disp([[0:length(nAD)-1]' rW rM nAD]);
% Write to excel to export to word
table1=[[0:length(nAD)-1]' rW rM nAD]; xlswrite('table2.xls',table1);