% Using data from the 'SVGS Specs June 2020' document, sets variance
% deviations for the error in all lateral and rotational degrees of freedom
% based on the true distance between the target and sensor. 

% R is defined for the state vector [X, Y, Z, Roll, Pitch, Yaw]
function [R] = svgs_R_from_range_SRT(range)

% Assuming linear error with a slope of: 
% [x y z phi theta psi]
% x = [0.0515; 0.0515; 0.018; 0.1324; 0.1324; 0.1324]; % Degrees
x = [0.0515; 0.0515; 0.018; 0.0023; 0.0023; 0.0023]; % Radians
% x = [0.0075; 0.0075; 0.0075; 0.0075; 0.0075; 0.0075]; % 5% of distance

% Slope values are for 3-sigma error, so dividing by 3
Range = (range / 3) * eye(6);
R_std = Range.*x;
% Compute variance from standard deviation
R = R_std*R_std;

end