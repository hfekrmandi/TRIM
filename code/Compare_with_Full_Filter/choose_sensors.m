function [C,ind_sensor] = choose_sensors(Ct,nSensor)
% This function chooses nSensor rows from Ct that have the highest sum
% of absolute value. The logic is that we want to choose the lines that 
% provide meaningful observation matrix.
norm_c = zeros(size(Ct,1),1);
for i = 1 : size(Ct,1)
    norm_c(i) = sum(abs(Ct(i,:)));
end

[B,ind] = sort(norm_c,'descend');
ind_sensor = ind(1:nSensor);
C = Ct(ind(1:nSensor),:);