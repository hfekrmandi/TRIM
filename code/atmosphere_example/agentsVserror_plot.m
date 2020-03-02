clearvars
close all
clc
load('/home/naveed/Documents/DSE_data/10to50.mat')
agents_count = [10 20 30 40 50];

T = size(error_(10).e_vs_gold,2);    
N = size(agents_count,2);
Hyb_errors = zeros(T, 2*N);
ICI_errors = zeros(T, 2*N);
for i=1:N
    Hyb_errors(:,i) = error_(agents_count(i)).e_vs_gold(1,:);
    ICI_errors(:,i) = error_(agents_count(i)).e_vs_gold(2,:);
end

load('/home/naveed/Documents/DSE_data/60to100.mat')
agents_count = [60 70 80 90 100];

for i=1:N
    Hyb_errors(:,i+N) = error_(agents_count(i)).e_vs_gold(1,:);
    ICI_errors(:,i+N) = error_(agents_count(i)).e_vs_gold(2,:);
end
agents_count = [10 20 30 40 50 60 70 80 90 100];
boxplot(Hyb_errors, agents_count)%, 'color','r')
hLegend = legend(findall(gca,'Tag','Box'), {'Hybrid'});
hold on 
boxplot(ICI_errors, agents_count,'PlotStyle','compact')
hLegend = legend(findall(gca,'Tag','Box'), {'ICI'});
ylim([0 0.06])