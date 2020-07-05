clear all

rng('default')

n_states = 2;
n_agents = 3;
n_obs = 1;
n_oa = n_obs*factorial(n_agents);
n_sa = n_states*n_agents;

x_prior = (1:n_states*n_agents)';
S_prior = 10*eye(n_sa);
y_prior = S_prior*x_prior;

x_11 = x_prior;
P_11 = inv(S_prior);
Y_11 = S_prior;
y_11 = S_prior*x_prior;

Q_i = 0.01*eye(n_states);
R_i = 0.01*eye(n_obs);
x_true = [3; 0; -1; 3; 5; 5];
F = eye(n_sa);
Q = 0.01*eye(n_sa);
R = 0.01*eye(n_oa);

% Add convergence plot

for i = 1:100
    x_true = F*x_true + mvnrnd(zeros(1, n_sa), Q)';
%     z_val = sqrt((x_true(3) - x_true(1))^2 + (x_true(4) - x_true(2))^2);
%     z = [z_val; z_val] + mvnrnd(zeros(1, n_oa), R)';
%     
%     H = [(x_11(1) - x_11(3))/z(1), (x_11(2) - x_11(4))/z(1), ...
%          (x_11(3) - x_11(1))/z(1), (x_11(4) - x_11(2))/z(1);...
%          (x_11(1) - x_11(3))/z(2), (x_11(2) - x_11(4))/z(2),...
%          (x_11(3) - x_11(1))/z(2), (x_11(4) - x_11(2))/z(2)];
    
    z_test = [];
    H_test = [];
    for i = 1:n_agents
        for j = [1:i-1, i+1:n_agents]
            z_temp = z_range_2d(x_true, i, j);
            z_test = vertcat(z_temp, z_test);
            
            H_temp = Hz_range_2d(x_11, i, j);
            H_test = vertcat(H_temp, H_test);
        end
    end
    z = z_test + mvnrnd(zeros(1, n_oa), R)';
    H_test = H_test ./ z;
    H = H_test;
    SBZ = z - H*x_11;
    
    M = inv(F)'*Y_11*inv(F);
    C = M*inv(M+inv(Q));
    L = eye(n_sa) - C;
    Y_01 = L*M*L' + C*inv(Q)*C';
    y_01 = L*inv(F)'*y_11;
    x_01_inf = inv(Y_01)*y_01;
    P_01_inf = inv(Y_01);
    Y_00 = Y_01 + H'*inv(R)*H;
    y_00 = y_01 + H'*inv(R)*z;
    x_00_inf = inv(Y_00)*y_00;
    P_00_inf = inv(Y_00);
    
    x_01 = F*x_11;
    P_01 = F*P_11*F' + Q;
    y = z - H*x_01;
    S = H*P_01*H' + R;
    K = P_01*H'*inv(S);
    x_00 = x_01 + K*y;
    P_00 = (eye(n_sa) - K*H)*P_01;
    y = z - H*x_00;
    
    x_11 = x_00;
    P_11 = P_00;
    
    Y_11 = Y_00;
    y_11 = y_00;
    
    delta_x_01 = x_01 - x_01_inf;
    delta_P_01 = P_01 - P_01_inf;
    delta_x_00 = x_00 - x_00_inf;
    delta_P_00 = P_00 - P_00_inf;
end

% for i = 1:100
%     x_true = F*x_true + mvnrnd([0, 0, 0, 0], Q)';
%     z_val = sqrt((x_true(3) - x_true(1))^2 + (x_true(4) - x_true(2))^2);
%     z = [z_val; z_val] + mvnrnd([0, 0], R)';
%     H = [(x_prior(1) - x_prior(3))/z(1), (x_prior(2) - x_prior(4))/z(1), ...
%         (x_prior(3) - x_prior(1))/z(1), (x_prior(4) - x_prior(2))/z(1);...
%         (x_prior(1) - x_prior(3))/z(2), (x_prior(2) - x_prior(4))/z(2),...
%         (x_prior(3) - x_prior(1))/z(2), (x_prior(4) - x_prior(2))/z(2)];
%     SBZ = z - H*x_01;
% 
%     S = S_prior + H'*inv(R)*H;
%     y = y_prior + H'*inv(R)*z;
%     
%     S_post = inv(Q + F*inv(S)*F');
%     y_post = S_post*F*inv(S)*y;
%     
%     S_prior = S_post;
%     y_prior = y_post;
%     x_post = inv(S_post)*y_post;
%     x_prior = x_post;
% end

function [Hz] = Hz_range_2d(x, agent1, agent2)
    n = 2;

    agent1_x = n*(agent1 - 1) + 1;
    agent1_y = agent1_x + (n - 1);
    agent2_x = n*(agent2 - 1) + 1;
    agent2_y = agent2_x + (n - 1);

    n_states = size(x,1);
    set_1 = zeros(n_states);
    set_2 = zeros(n_states);

    set_1(agent1_x:agent1_y,agent1_x:agent1_y) = eye(n);
    set_1(agent2_x:agent2_y,agent2_x:agent2_y) = eye(n);
    set_2(agent1_x:agent1_y,agent2_x:agent2_y) = eye(n);
    set_2(agent2_x:agent2_y,agent1_x:agent1_y) = eye(n);

    Hz = (set_1*x - set_2*x)';
end

function [z] = z_range_2d(x, agent1, agent2)
    n = 2;

    agent1_x = n*(agent1 - 1) + 1;
    agent1_y = agent1_x + (n - 1);
    agent2_x = n*(agent2 - 1) + 1;
    agent2_y = agent2_x + (n - 1);

    n_states = size(x,1);
    set_1 = zeros(n_states);
    set_2 = zeros(n_states);

    set_1(agent1_x:agent1_y,agent1_x:agent1_y) = eye(n);
    %set_1(agent2_x:agent2_y,agent2_x:agent2_y) = eye(n);
    set_2(agent1_x:agent1_y,agent2_x:agent2_y) = eye(n);
    %set_2(agent2_x:agent2_y,agent1_x:agent1_y) = eye(n);
    
    z = sqrt((set_1*x - set_2*x)'*(set_1*x - set_2*x));
end