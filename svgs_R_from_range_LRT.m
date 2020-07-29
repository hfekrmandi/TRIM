% Using data from the 'SVGS Specs June 2020' document, sets standard
% deviations for the error in all lateral and rotational degrees of freedom
% based on the true distance between the target and sensor. 

% R is defined for the state vector [X, Y, Z, Roll, Pitch, Yaw]
function [R] = svgs_R_from_range_LRT(range)

% Assuming linear error with a slope of: 
x = [0.0666; 0.001; 0.001; 0.0074; 0.0349; 0.0349];

% Slope values are for 3-sigma error, so dividing by 3
Range = (range / 3) * eye(6);
R_std = Range.*x;
% Compute variance from standard deviation
R = R_std*R_std;

end