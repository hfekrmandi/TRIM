% Run the hypothesizing distributed Kalman filter in a simple scenario
% @param NUM_SENSORS			number of sensors in the network
% @param NUM_TIMESTEPS			number of time steps the model is evaluated
% @param PACKET_DROP_RATE		probability that transmission from sensor to fusion node fails
% @param SENSOR_RESTART_PROB	probability that sensor has malfunction and is restarted
% @param K_UPDATE_FREQUENCY		number of time steps after the HGMM is adaptively updated (set to -1 to disable)
% @param K_UPDATE_DELAY			number of time steps until a new calculated HGMM is employed.
function RunExample(NUM_SENSORS, NUM_TIMESTEPS, PACKET_DROP_RATE, SENSOR_RESTART_PROB, K_UPDATE_FREQUENCY, K_UPDATE_DELAY)
	rng('shuffle');
	
	SENSOR_RESTART_PROB = SENSOR_RESTART_PROB / 100;
	PACKET_DROP_RATE = PACKET_DROP_RATE / 100;
		
	% suppress spam messages
	if(SENSOR_RESTART_PROB > 0)
		warning('off','HKF:CX_CHECKSUM')
	end

	xTrue = [20;20];
	A = [1 1 ; 0 1];
	Cw = 10 * eye(2);
	H = eye(2);
	Cv = eye(2);
	CxInit = eye(2);
	CzInit = eye(2);

	%
	% Declaration
	%
	network = SimpleNetwork;
	fusion_node = HKF_FusionNode;
	sensors = cell(NUM_SENSORS);

	%
	% Initialization
	%
	fusion_node.init(network, 0, NUM_SENSORS, CxInit, CzInit, 10, K_UPDATE_FREQUENCY, K_UPDATE_DELAY);
	for i=1:NUM_SENSORS
		sensors{i} = HKF_Sensor;
		sensors{i}.init(network, i, 0, CxInit, CzInit);
	end
	network.init(sensors, fusion_node);

	%
	% Simple estimation scenario
	%
	xErr = zeros(NUM_TIMESTEPS,1);
	cxApprox = zeros(NUM_TIMESTEPS,1);
	cxBound = zeros(NUM_TIMESTEPS,1);
	for k = 1 : NUM_TIMESTEPS
		xTrue = A * xTrue + mvnrnd(zeros(size(Cw,1),1),Cw)';
		
		% process sensors and fusion node
		fusion_node.process(A, Cw);
		for i=1:NUM_SENSORS
			sensors{i}.process(A, Cw);
		end    
		
		% simulate measurements    
		for i=1:NUM_SENSORS
			z = H * xTrue + mvnrnd(zeros(size(H,1),1),Cv)';
			sensors{i}.addMeasurement(z, H, Cv);
		end
		
		% simulate communication
		for i=1:NUM_SENSORS
			if rand() > 1 - PACKET_DROP_RATE
				continue;
			end
			sensors{i}.sendNetworkPacket();
		end

		% simulate sensor crashes
		for i=1:NUM_SENSORS
			if rand() > SENSOR_RESTART_PROB
				continue;
			end
			% technique: init with wrong CxInit -> restart is triggerd by
			% fusion node
			sensors{i}.init(network, i, k, CxInit, CzInit);
		end
		
		% set errors
		[x, CxApprox, CxBound] = fusion_node.getEstimate();
		xErr(k) = (xTrue - x)'*(xTrue - x);
		cxApprox(k) = trace(CxApprox);
		cxBound(k) = trace(CxBound);
	end

	%
	% Evaluation
	%
	plot(1:NUM_TIMESTEPS, xErr, ...
		1:NUM_TIMESTEPS, cxApprox,...
		1:NUM_TIMESTEPS, cxBound);
end