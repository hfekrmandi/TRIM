clearvars
close all
clc
load('/home/naveed/Documents/DSE_data/80_states_convergence_rate_30_agents.mat')
converg_steps_array = [1,5,10,15,20,30,40,50,60,80,100,120];

for i=1:9
    error_new(i) = error_(i);
end

load('/home/naveed/Documents/DSE_data/80_states_convergence_rate_30_agents_80to120.mat')
for i=10:12
    error_new(i) = error_(i-9);
end
converg_steps_array = [1,5,10,15,20,30,40,50,60,80,100,120];

clearvars -except error_new converg_steps_array