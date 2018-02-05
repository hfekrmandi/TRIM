classdef HKF_Sensor < handle
	properties (Access = protected)
        % network
		g_network_object;	% an implementation of HKF_Network
        
		% common
		sensor_idx;		% index of the sensor
		k_cur;			% current time step
		Cx;				% auxiliary variable ("global covariance")
		Cz;				% HGMM
		x;				% sensor estimate
		Delta;			% delta matrix
		initalized;		% is set to one when at least one measurement or estimate has been added
		% % %
		% OPTIONAL: bound calculations
		Bx;				% Bx matrix
		Bw;				% Bw matrix
		
		% % %
		% OPTIONAL: adaptive improvement of HGMM
		k_next_update;	% time step when HGMM is replaced. -1 for never
		Cz_New;			% new HGMM
	end
	
    methods	
		%
		% Initialize the sensor node
		% 
		% @param network_object 	an implementation of HKF_Network
        % @param sensor_idx unique index of the sensor node
		% @param k_cur		current time step
		% @param CxInit		an arbitrary covariance matrix that must initially be the same at all nodes
		%					and should reflect the steady state covariance of the KF with global
		%					measurement model CzInit best possible.
		% @param CzInit		initial value for the HGMM
		function init(o, network_object, sensor_idx, k_cur, CxInit, CzInit)
			dimX = length(CxInit);
			
            o.g_network_object = network_object;
            
			o.initalized = 0;
			o.sensor_idx = sensor_idx;
			o.k_cur = k_cur;
			o.Cx = CxInit;
			o.Cz = CzInit;
			
			o.x = zeros(dimX,1);
			o.Delta = zeros(dimX,dimX);
			o.Bx = zeros(dimX,dimX);
			o.Bw = zeros(dimX,dimX);
			
			o.k_next_update = -1;
		end
	
		%
		% Predict the available estimate to the next time step
		%
		% @param A			state transition matrix
		% @param Cw			state transition uncertainty
		function process(o, A, Cw)
            dimX = length(Cw);
			o.k_cur = o.k_cur + 1;
			
			% update cz
			if(o.k_next_update~=-1 && o.k_cur == o.k_next_update)
				o.Cz = o.Cz_New;
			end
			
			CxPredInv = eye(dimX)/(A * o.Cx * A' + Cw);
			o.Cx = eye(dimX)/(CxPredInv + o.Cz);
			M = o.Cx * CxPredInv * A;
			A_Inv = eye(dimX)/A;
			
			% update variables
			o.x = M * o.x ; 
			o.Delta = M * o.Delta * A_Inv;
            o.Bx = M * o.Bx * M';
			o.Bw = M * o.Bw * M' + o.Delta * Cw * o.Delta';
		end

	
		%
		% Add estimate information to the local estimate
		% 
		% @param xInit 		initial estimate
		% @param CxEstInit 	initial covariance matrix
		function addEstimate(o, xInit, CxEstInit)   
			dimX = length(CxEstInit);
			o.addMeasurement(xInit, eye(dimX), CxEstInit);
		end
		
		%
		% Add the information of a new measurement to the local estimate
		%
		% @param z 			measurement
		% @param H 			measurement matrix
		% @param Cv			measurement noise covariance matrix
		function addMeasurement(o, z, H, Cv)
			L = o.Cx * H' / Cv;
			
			o.initalized = 1;
			o.x = o.x + L * z;
			o.Delta = o.Delta + L * H;
			o.Bx = o.Bx  + L * Cv * L';
		end
		
		%
		% Update the HGMM
		%
		% @param k_to_update 	time step when new HGMM should be used the first time
		% @param CzNew 			new HGMM
		function updateHGMM(o, k_to_update, CzNew)
			o.k_next_update = k_to_update;
			o.Cz_New		= CzNew;
		end
		
		%
		% Send data to the fusion node
		%
		function sendNetworkPacket(o)
			if o.initalized == 0
				warning('HKF:NOT_INITIALIZED','Could not send sensor estimate. No information contained.');
				return;
			end
            
            o.g_network_object.sendToFusionNode(o.sensor_idx, o.k_cur, o.x, ...
                o.Delta, HKF_Network.getchecksum(o.Cx), o.Bx, o.Bw);
			end
	end
end