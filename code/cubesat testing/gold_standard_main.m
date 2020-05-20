%%main script
clear all
save_loc = 'cubesat_test_3.mat';
gold_standard();
save(save_loc);
%addpath('/libs/latexfigure/');
agents_state_plot(save_loc);