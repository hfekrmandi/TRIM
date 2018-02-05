classdef HKF_FusionNode < handle
	properties (Access = protected)
		NUM_TS_TO_STORE    = 10;  % max time steps old models are stored for out-of-sequence estimates
		K_UPDATE_FREQUENCY = 10;  % set to -1 to prevent updating
		K_UPDATE_DELAY	   = 3;   % time steps between calculation and application of new HGMM

		% network
		g_network_object;	% an implementation of HKF_Network

		% common
		k_cur;			% current time step
		Cx;				% auxiliary variable ("global covariance")
		Cz;				% HGMM
		x_vec;			% sensor estimates
		Delta_vec;		% delta matrices of sensors
		k_vec;			% time step when estimate has been sent by sensor
		init_vec;		% helper variable to store which estimates are initialized
		
		% % %
		% OPTIONAL: bound calculations
		Bx_vec;		% Bx matrices of sensors
		Bw_vec;		% Bw matrices of sensors

		% % %
		% OPTIONAL: out-of-sequence estimates
		% matrices of time steps k_hist_begin to k_cur-1 are stored in matrices 
		k_hist_begin;	 % models are available for k_hist_begin+1
		cx_checksum_hist;% Cx checksum history 
		M_hist;			 % M = K_{t+1} * A_t 
		Cw_hist			 % Cw_t 
		A_Inv_hist;		 % inv{A_t} 
		
		% % %
		% OPTIONAL: adaptive improvement of HGMM
		Cz_New;
	end
    
    methods
		%
		% Initialize the fusion node
		% 
		% @param network_object 	an implementation of HKF_Network
		% @param k_cur          current time step
		% @param nNodes         number of sensors in the network
		% @param CxInit         an arbitrary covariance matrix that must initially be the same at all nodes
		%                       and should reflect the steady state covariance of the KF with global
		%                       measurement model CzInit best possible.
		% @param CzInit         initial value for the HGMM
        % @param varargin{1}    number of time steps old models are stored
        % @param varargin{2}    number of time steps after the HGMM is
        %                       adaptively updated (set to -1 to disable)
        % @param varargin{3}    number of time steps until a new calculated
        %                       HGMM is employed.
        function init(o, network_object, k_cur, nNodes, CxInit, CzInit, varargin)   
			dimX = length(CxInit);
			
			o.g_network_object = network_object;
			o.k_cur = k_cur;
			o.Cx = CxInit;
			o.Cz = CzInit;
			o.x_vec = zeros(dimX*nNodes,1);
			o.Delta_vec = zeros(dimX*nNodes,dimX);
			o.k_vec = -ones(nNodes,1);
			o.init_vec = zeros(nNodes, 1);
			
			o.Bx_vec = zeros(dimX*nNodes,dimX);
			o.Bw_vec = zeros(dimX*nNodes,dimX);
			
			o.k_hist_begin = k_cur;
			o.cx_checksum_hist = {}; 
			o.M_hist = {}; 
			o.Cw_hist = {}; 
			o.A_Inv_hist = {}; 
            
			o.Cz_New = CzInit;
            
            % set optional parameters
            if(~isempty(varargin)) 
                o.NUM_TS_TO_STORE = varargin{1}; end
            if(length(varargin)>1) 
                o.K_UPDATE_FREQUENCY = varargin{2}; end
            if(length(varargin)>2) 
                o.K_UPDATE_DELAY = varargin{3}; end
        end
        
		%
		% Add an estimate of a sensor to the fusion node
		%
		% @param sensor_idx		sensor id
		% @param k 				time stamp of estimate
		% @param x 				estimate
		% @param Delta 			correction matrix
		% @param varargin{1} 	cx checksum
		% @param varargin{2} 	Bx
		% @param varargin{3} 	Bw
		function addEstimate(o, sensor_idx, k, x, Delta, varargin)
			% check if estimate is older than an earlier received one
			if k <= o.k_vec(sensor_idx)
				warning('HKF:OLD_MSG','Received old measurement. Drop...');
				return;
			end
			
			% check if our out of sequence processing has enough data
            if k ~= o.k_cur && k <= o.k_hist_begin
				warning('HKF:OOS_LIMITED','Cannot include estimate. Old models have been dropped');
				return;
            elseif k > o.k_cur
                warning('HKF:OOS_LIMITED','Cannot include estimate. Fusion node time smaller than estimate time');
                return;
            end
			
            
			% check if checksum code of estimate matches the one of the fusion node
			cx_checksum = HKF_Network.getchecksum(o.Cx);
            if(k < o.k_cur)
                cx_checksum = o.cx_checksum_hist{k - o.k_hist_begin + 1};
            end
			if ~isempty(varargin) && abs(varargin{1}-cx_checksum) > 1e-6
				warning('HKF:CX_checksum','Cx checksum does not match with fusion center. k = %d',k);
				o.g_network_object.restartSensor(sensor_idx, o.k_cur, o.Cx, o.Cz);
				return;
			end				
			
			% set estimate state to initialized
			o.init_vec(sensor_idx) = 1;
			
			% handle optional parameters
			dimX = length(Delta);
			Bx = zeros(dimX);
			Bw = zeros(dimX);
			if(length(varargin)>1)
				Bx = varargin{2};
			end 
			if(length(varargin)>2)
				Bw = varargin{3};
			end

			% predict to current time step
			[x, Delta, Bx, Bw] = o.predict(k, o.k_cur, x, Delta, Bx, Bw);
			 
			% update parameters			
			xIdx = dimX * (sensor_idx-1) + 1 : dimX * sensor_idx;
			o.x_vec(xIdx) = x; 
			o.Delta_vec(xIdx,:) = Delta;
			o.k_vec(sensor_idx) = k;
			o.Bx_vec(xIdx,:) = Bx;
			o.Bw_vec(xIdx,:) = Bw;
		end
		
		%
		% Predict the available estimates to the next time step and
		% calculate a new hypothesis when necessary
		%
		% @param A				state transition matrix
		% @param Cw				state transition uncertainty
        function process(o, A, Cw)
			% predict + filter
			dimX = length(Cw);
			nEstimates = length(o.x_vec)/dimX;
                       
            % set new hypothesis (based on filtered (not predicted) data)
			if(o.K_UPDATE_FREQUENCY~=-1 && mod(o.k_cur, o.K_UPDATE_FREQUENCY) == 0 && sum(o.init_vec) > 0)
                DeltaFus = repmat(eye(dimX),1,nEstimates) * o.Delta_vec;
                deltaFusVector = reshape(DeltaFus, dimX*dimX, 1);
                M = (eye(dimX)-o.Cx*o.Cz)*A;
                Gain = eye(dimX*dimX) - kron(eye(dimX)/A', M);
                o.Cz_New = eye(dimX)/o.Cx * reshape(Gain * deltaFusVector, dimX, dimX);
				o.g_network_object.broadcastNewHGMM(o.k_cur + o.K_UPDATE_DELAY, o.Cz_New);
			end
            
            % update hypothesis when it is time
			o.k_cur = o.k_cur + 1;
            if(o.K_UPDATE_FREQUENCY ~= -1 && ...
                    mod(o.k_cur,o.K_UPDATE_FREQUENCY) == o.K_UPDATE_DELAY)
				o.Cz = o.Cz_New;
            end
			
            CxOld = o.Cx;
            CxPred = A * o.Cx * A' + Cw;
			o.Cx = eye(dimX)/(eye(dimX)/CxPred + o.Cz);
			M = o.Cx / CxPred * A; % K * A;
			A_Inv = eye(dimX)/A;
            
			% set history variables
			if(length(o.cx_checksum_hist) >= o.NUM_TS_TO_STORE)
				o.k_hist_begin = o.k_hist_begin + 1;
				o.cx_checksum_hist(1) = [];
				o.M_hist(1) = [];
				o.Cw_hist(1) = [];
				o.A_Inv_hist(1) = [];
			end
			newIdx = length(o.cx_checksum_hist) + 1;
			o.cx_checksum_hist{newIdx} = HKF_Network.getchecksum(CxOld);
			o.M_hist{newIdx} = M;
			o.Cw_hist{newIdx} = Cw;
			o.A_Inv_hist{newIdx} = A_Inv;
	
			% predict estimates to current time step
            for i = 1 : nEstimates
                if o.init_vec(i) == 0 % not yet initialized
					continue;
                end
				xIdx = dimX * (i-1) + 1 : dimX * i;
				[o.x_vec(xIdx), o.Delta_vec(xIdx,:), o.Bx_vec(xIdx,:), o.Bw_vec(xIdx,:)] = o.predict(o.k_cur-1, ...
					o.k_cur, o.x_vec(xIdx), o.Delta_vec(xIdx,:), o.Bx_vec(xIdx,:), o.Bw_vec(xIdx,:));
            end
		end   
	   
	    %
		% Return the fused estimate
	    function [x, CxApprox, CxBound] = getEstimate(o)
			dimX = length(o.Cx);
			nEstimates = length(o.x_vec)/dimX;
			
			if sum(o.init_vec) == 0
				error('HKF:NOT_INITIALIZED','No estimate has been initialized. Cannot obtain fused estimate');
			end
			
			x_fus = zeros(dimX,1);
			Delta_fus = zeros(dimX,dimX);
			Bx_fus = zeros(dimX,dimX);
			Bw_fus = zeros(dimX,dimX);
			for i=1:nEstimates
				if(o.init_vec(i) == 0) % not yet initialized
					continue;
				end
				xIdx = dimX * (i-1)+1:dimX * i;
				x_fus = x_fus + o.x_vec(xIdx);
				Delta_fus = Delta_fus + o.Delta_vec(xIdx,:);
				Bx_fus = Bx_fus + o.Bx_vec(xIdx,:);
				Bw_fus = Bw_fus + o.Bw_vec(xIdx,:);
			end
			
			if(cond(Delta_fus)> 1e+15)
				error('HKF:NOT_INITIALIZED','Delta matrix is close to singular. Cannot obtain estimate.');
			end
			Delta_fus_Inv = eye(dimX) / Delta_fus;
			
			% Bw approximation
			tLastIdx = length(o.cx_checksum_hist);
			M = o.M_hist{tLastIdx};
			N = M * Delta_fus * o.A_Inv_hist{tLastIdx} * o.Cw_hist{tLastIdx} * ...
				o.A_Inv_hist{tLastIdx}' * Delta_fus' * M';
			normM = max(max(M));
			Bw_approx = N;
			for i = 1 : o.k_cur - 1
				normRest = max(max(N)) * normM.^2 * (1 - normM.^(2*o.k_cur - 2*i)) / (1 - normM.^2) ;
				normBw = max(max(Bw_approx));
				if(normRest / normBw   < 0.005)
					break;
				end
				N = M * N * M';				
				Bw_approx = Bw_approx + N;
			end
			
			% set results
			x = Delta_fus_Inv * x_fus;
			CxApprox = Delta_fus_Inv * (Bx_fus + Bw_approx) * Delta_fus_Inv';
			CxBound = Delta_fus_Inv * (Bx_fus + nEstimates * Bw_fus) * Delta_fus_Inv';
		end
	end

	methods(Access = protected)
       function [x, Delta, Bx, Bw] = predict(o, k_from, k_to, x, Delta, Bx, Bw)
           if(k_from == k_to)
               return;
           end
           
           if k_from < o.k_hist_begin
               error('HKF:OOS_LIMITED','Could not predict estimates. OOS data not old enough.');
            end
            
            tIdx = k_from - o.k_hist_begin + 1;
            M = o.M_hist{tIdx};
            Cw = o.Cw_hist{tIdx};
            A_Inv = o.A_Inv_hist{tIdx};
			for k_temp = k_from : k_to - 1
				x = M*x;
				Delta = M * Delta * A_Inv ;
				Bx = M* Bx *M';
                Bw = M * Bw * M' + Delta * Cw * Delta';
			end
       end
	end
end