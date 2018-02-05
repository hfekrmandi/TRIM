function [x_est, P_est] = loc_kalman_update(x_pred,P_pred,z,h,H)
global opt_dist

% Update aposteriori state estimate
R = opt_dist.obs.R*eye(size(H,1));
K = P_pred*H'*(H*P_pred*H'+R)^(-1);
x_est = x_pred + K * (z - h);
% z - h;
% Update aposteriori error covariance estimate
P_est = (eye(length(x_pred)) - K*H) * P_pred;


end