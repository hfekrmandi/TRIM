
OMAS_initialize
{
	% Set up simulation data and define agents and their public values
	ConfigureSimulationParameters
	% Simulate the system
	OMAS_process
	{
		% Simulate time steps
		ProcessGlobalTimeSeries
		{
			for each time step:
			{
				if agents are not idle:
				{
					% Store global states
					% Compute states at the next time step
					OMAS_updateGlobalStates

					Step 6: COMPUTE INFORMATION FILTER ESTIMATION (@t=k) - OMAS_process Line 124
					UpdateInformationFilter
					{
						% Measurement = relative (x, y) position + noise
						observedz = observedPosition + mvnrnd(zeros(1, 2), noise);

						% Pass values to the agents information filter
						referenceObject.InformationFilter - agent.m Line 240
						{
							% Compute information filter steps
							% Store measurements
							% Store motion model predictions
						}
					}

					Step 7: COMPUTE CONSENSUS STEPS (@t=k) - OMAS_process Line 147
					apply_comm_model - OMAS_process Line 149
					{
						% Determine which agents can communicate with one another
					}
					% break agents into groups of agents that can communicate with eachother
					consensus - OMAS_process Line 151
					{
						for step = 1:20
						{
							consensus_group - OMAS_process Line 800
							{
								% Compute consensus of motion model predictions - OMAS_process Line 825
								% Compute consensus of measurements - OMAS_process Line 830
								% Compute consensus - OMAS_process Line 838 & 839
							}
						}
					}

					Step 8: UPDATE AGENT ESTIMATE FROM CONSENSUS DATA (@t=k) - OMAS_process Line 154
					UpdateObject
					{
						% Measurement = consensus estimate

						% Pass values to the agents information filter
						referenceObject.main - agent.m Line 138
						{
							% Compute control signal
						}
					}
				}
			}
		}
	}
}