function plot_e_est(e_est,idx_agent)
global opt_dist

plot(squeeze(e_est((idx_agent-1)*(opt_dist.dimAgents)+1,1,:)),'b');hold on
plot(squeeze(e_est((idx_agent-1)*(opt_dist.dimAgents)+1,2,:)),'r')
plot(squeeze(e_est((idx_agent-1)*(opt_dist.dimAgents)+1,3,:)),'g')
plot(squeeze(e_est((idx_agent-1)*(opt_dist.dimAgents)+1,4,:)),'k')
legend(['1\rightarrow e_{',num2str(idx_agent),'}'],...
    ['2\rightarrow e_{',num2str(idx_agent),'}'],...
    ['3\rightarrow e_{',num2str(idx_agent),'}'],...
    ['4\rightarrow e_{',num2str(idx_agent),'}']);
plot(squeeze(e_est((idx_agent)*(opt_dist.dimAgents),1,:)),'--b');hold on
plot(squeeze(e_est((idx_agent)*(opt_dist.dimAgents),2,:)),'--r')
plot(squeeze(e_est((idx_agent)*(opt_dist.dimAgents),3,:)),'--g')
plot(squeeze(e_est((idx_agent)*(opt_dist.dimAgents),4,:)),'--k')
end