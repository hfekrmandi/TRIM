function plot_x_est(x_est,idx_agent)
global opt_dist

plot(squeeze(x_est((idx_agent-1)*(opt_dist.dimAgents)+1,1,:)),'b');hold on
plot(squeeze(x_est((idx_agent-1)*(opt_dist.dimAgents)+1,2,:)),'r')
plot(squeeze(x_est((idx_agent-1)*(opt_dist.dimAgents)+1,3,:)),'g')
plot(squeeze(x_est((idx_agent-1)*(opt_dist.dimAgents)+1,4,:)),'k')
legend(['1\rightarrow x_{',num2str(idx_agent),'}'],...
    ['2\rightarrow x_{',num2str(idx_agent),'}'],...
    ['3\rightarrow x_{',num2str(idx_agent),'}'],...
    ['4\rightarrow x_{',num2str(idx_agent),'}']);
plot(squeeze(x_est((idx_agent)*(opt_dist.dimAgents),1,:)),'--b');hold on
plot(squeeze(x_est((idx_agent)*(opt_dist.dimAgents),2,:)),'--r')
plot(squeeze(x_est((idx_agent)*(opt_dist.dimAgents),3,:)),'--g')
plot(squeeze(x_est((idx_agent)*(opt_dist.dimAgents),4,:)),'--k')

end