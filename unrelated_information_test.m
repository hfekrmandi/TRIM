rng('default')

x_prior = [0; 1; 2; 3];
S_prior = 2*eye(4);
y_prior = S_prior*x_prior;

Y_11 = S_prior;
y_11 = S_prior*x_prior;

Q_i = 0.1*eye(2);
R_i = 0.1;
x_true = [3; 0; -1; 3];
F = eye(4);
Q = 0.1*eye(4);
R = 0.1*eye(2);

for i = 1:100
    x_true = F*x_true + mvnrnd([0, 0, 0, 0], Q)';
    z_val = sqrt((x_true(3) - x_true(1))^2 + (x_true(4) - x_true(2))^2);
    z = [z_val; z_val] + mvnrnd([0, 0], R)';
    
    x_11 = inv(Y_11)*y_11;
    
    H = [(x_11(1) - x_11(3))/z(1), (x_11(2) - x_11(4))/z(1), ...
         (x_11(3) - x_11(1))/z(1), (x_11(4) - x_11(2))/z(1);...
         (x_11(1) - x_11(3))/z(2), (x_11(2) - x_11(4))/z(2),...
         (x_11(3) - x_11(1))/z(2), (x_11(4) - x_11(2))/z(2)];
    SBZ = z - H*x_11;

    M = inv(F)'*Y_11*inv(F);
    C = M*inv(M+inv(Q));
    L = eye(4) - C;
    Y_01 = L*M*L' + C*inv(Q)*C';
    y_01 = L*inv(F)'*y_11;
    Y_00 = Y_01 + H'*inv(R)*H;
    y_00 = y_01 + H'*inv(R)*z;
    
    Y_11 = Y_00;
    y_11 = y_00;
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
