classdef SimpleNetwork < HKF_Network 
	properties (Access = protected)
		fusion_node;	% fusion node
		sensors;		% cell array with sensors
	end
	
	methods
		function init(obj, sensors, fusion_node)
			obj.sensors = sensors;
			obj.fusion_node = fusion_node;
		end
	
		function sendToFusionNode(obj, sensor_idx, k_cur, x, Delta, varargin)
			obj.fusion_node.addEstimate(sensor_idx, k_cur, x, Delta, varargin{:});
		end
		
		function restartSensor(obj, sensor_idx, k_cur, CxInit, Cz)
			obj.sensors{sensor_idx}.init(obj, sensor_idx, k_cur, CxInit, Cz);
		end
		
		function broadcastNewHGMM(obj, k_update, CzNew)
			for i=1:length(obj.sensors)
				obj.sensors{i}.updateHGMM(k_update, CzNew);
			end
		end
	end
end