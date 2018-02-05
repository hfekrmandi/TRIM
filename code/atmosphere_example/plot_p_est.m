function plot_p_est(P_est,idx_agent)
% global opt_dist

plot(sqrt(squeeze(P_est(1,1,idx_agent,:))),'b');hold on
plot(sqrt(squeeze(P_est(3,3,idx_agent,:))),'r')
plot(sqrt(squeeze(P_est(5,5,idx_agent,:))),'g')
plot(sqrt(squeeze(P_est(7,7,idx_agent,:))),'k')
legend([num2str(idx_agent),'\rightarrow\sigma_{11}'],...
    [num2str(idx_agent),'\rightarrow\sigma_{22}'],...
    [num2str(idx_agent),'\rightarrow\sigma_{33}'],...
    [num2str(idx_agent),'\rightarrow\sigma_{44}']);

plot(sqrt(squeeze(P_est(2,2,idx_agent,:))),'.-b');hold on
plot(sqrt(squeeze(P_est(4,4,idx_agent,:))),'.-r')
plot(sqrt(squeeze(P_est(6,6,idx_agent,:))),'.-g')
plot(sqrt(squeeze(P_est(8,8,idx_agent,:))),'.-k')


end