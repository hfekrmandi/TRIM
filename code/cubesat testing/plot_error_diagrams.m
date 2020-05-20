function plot_error_diagrams(dataG)
global opt_dist
i_time = opt_dist.i_time;
figure(dataG.h_error);
subplot(411)
for i_error=1:opt_dist.dimState
    plot(squeeze(dataG.e_est(i_error,9,2:i_time))) ; hold on;
    title(['est error 1'])
end
subplot(412)
for i_error=1:opt_dist.dimState
    plot(squeeze(dataG.e_est(i_error,10,2:i_time))) ; hold on;
    title(['est error 2'])
end
subplot(413)
for i_error=1:opt_dist.dimState
    plot(squeeze(dataG.e_est(i_error,11,2:i_time))) ; hold on;
    title(['est error 3'])
end
subplot(414)
for i_error=1:opt_dist.dimState
    plot(squeeze(dataG.e_est(i_error,16,2:i_time))) ; hold on;
    title(['est error 4'])
end



figure(dataG.h_error2);
subplot(411)
for i_error2=1:opt_dist.nAgents
    plot((dataG.e_i_con(i_error2,4:i_time))) ; hold on;
    title(['norm of i\_{con} error 1'])
end
subplot(412)
for i_error2=1:opt_dist.nAgents
    plot((dataG.e_I_con(i_error2,4:i_time))) ; hold on;
    title(['norm of I\_{con} error 1'])
end
subplot(413)
for i_error2=1:opt_dist.nAgents
    plot((dataG.e_est_norm(i_error2,4:i_time))) ; hold on;
    title(['norm of e\_{est} error 1'])
end

subplot(414)
for i_error2=1:opt_dist.nAgents
    plot((dataG.e_P_est_norm(i_error2,4:i_time))) ; hold on;
end
    title(['norm of e\_{cov} error 1'])



end