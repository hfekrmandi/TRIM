classdef HKF_Network < handle
	methods(Abstract)
		% cx_checksum, Bx, Bw
		sendToFusionNode(obj, sensor_idx, k_cur, x, Delta, varargin)
		restartSensor(obj, sensor_idx, k_cur, CxInit, Cz)
		broadcastNewHGMM(obj, k_update, CzNew)
    end
    
    methods (Static, Access = public)
        function checksum = getchecksum(M)
			checksum = size(M,1) * size(M,2) * sum(sum(M)) + min(max(M)) + max(min(M));
	    end 
    end
end