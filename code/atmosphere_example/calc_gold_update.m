function calc_gold_update()
global opt_dist
for i_agent = 1 : opt_dist.nAgents
    run_full_filter(i_agent)
end
end
function run_full_filter(i_agent)
global opt_dist



% initialization
x_bar = opt_dist.result.initial.x_bar(:,i_agent);
P_bar =  opt_dist.result.initial.P_bar(:,:,i_agent);

x_update = x_bar;
P_update = P_bar;

for i_step = 1 : opt_dist.i_step
    % predict
    [x_pred,P_pred] = f(x_update,P_update,0);
    Y_bar = pinv(P_pred);
    y_bar = Y_bar*x_pred;
    % calc update
    delta_I = zeros(size(Y_bar));
    delta_i = zeros(size(y_bar));
    
    [size_group,nComponents,members] = networkComponents_gold(opt_dist.Graph_History{i_step});
    neigbours=members{nComponents(i_agent)};
    
    
    for i_neighbour = 1:numel(neigbours)
        H = opt_dist.result.obs.H{i_step,neigbours(i_neighbour)};
        z = opt_dist.sim.obs.z{i_step,neigbours(i_neighbour)};
        delta_I = delta_I + 1/(opt_dist.sim.obs.r_var{i_step,neigbours(i_neighbour)})*(H'*H);
        delta_i = delta_i +  1/(opt_dist.sim.obs.r_var{i_step,neigbours(i_neighbour)})*H'*z;
    end
    % update
    Y_update = Y_bar + delta_I;
    y_update = y_bar + delta_i;
    
    P_update = pinv(Y_update);
    x_update = P_update*y_update;
end
opt_dist.result.est_gold{i_agent}.Y_bar = Y_update;
opt_dist.result.est_gold{i_agent}.y_bar = y_update;
opt_dist.result.est_gold{i_agent}.P_bar = P_update;
opt_dist.result.est_gold{i_agent}.x_bar = x_update;

x_update



end