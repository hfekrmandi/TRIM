% EX1 : Example 1. 
%       Script which calls opconstr.m to solve the digital optimal
%       control problem specified in dopcoex.m using automatic
%       differentiation and next computes an optimal LQG compensator

load ex1 % load data from optimal control example
hf=figure; plot(tt,xt(:,1:end-1)); hold on; stair(tu(:,1)',tu(1:end-1,2:end)','m');
hold off; xlabel('Time'); ylabel('States');
str=['Optimal control & state trajectory min(J)=' num2str(J)];
str=[str ' tf=' num2str(tu(end,1))]; title(str);
print(hf,'Fig1.emf','-dmeta'); % Save figure
%% Compute equivalent discrete-time LQG problem data
clear; lqgcomp('tlsdyn','tlsout','tlsphi','tlslqgc','dop1','montcar',3);
%% Temporal stabilizability and detectability analysis
% Compute equivalent discrete-time LQG problem data
% clear; edodata('tlsdyn','tlsout','tlsphi','tlslqgc','dop1')
load edorp; % Load equivalent discrete-time LQG problem data
js=2; ks=js; tolWM=1e-6; epsrw=1e-6; %tolWM=1e-9; epsrw=1e-9;
[npt,nst,ptt,stt,NW,NM,rW,rM,nAD,hfs]=tmpstdt(lqgarr,js,ks,tolWM,epsrw,epsrw);
print(hfs(1),'Fig2.emf','-dmeta'); % Save figure
print(hfs(2),'Fig3.emf','-dmeta'); % Save figure
% Display table
disp('     i rank(W) rank(M) nAD'); disp([[0:length(nAD)-1]' rW rM nAD]);
% Write to excel to export to word
table1=[[0:length(nAD)-1]' rW rM nAD]; xlswrite('table1.xls',table1);