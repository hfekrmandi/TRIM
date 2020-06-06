% Run the simulation
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

					% Update events and seperations
					UpdateSimulationMeta
					{
						for each agent:
						{
							for each other agent:
							{
								% compute relative position
								% check for collisions
							}
						}
						fdfdfdf
						% Check if each agent can see each other agent
							% Some form of geometry with centroid, normals, faces, and vertices
							% set detectionLogicals variable
						%% Need to add obstacles block LOS
						%% 
						for each agent:
						{
							for each other agent:
							{
								if detection:
								{
									% Add detection event
									OMAS_eventHandler
								}
								if waypoint reached:
								{
									% set waypoint to reached
									OMAS_eventHandler
								}
								if warning (close to obstacle):
								{
									% Add warning event
									OMAS_eventHandler
								}
								if collision:
								{
									% Add collision event
									OMAS_eventHandler
								}
							}
						}
					}
					for each object:
					{
						%
						UpdateObjects
						{
							for each detected object:
							{
								% create detection packet
								% Call main in agent with all packets
							}
						}
					}
				}
			}
		}
	}
}

Agents (In general)
Methods:
	main
	{
		Takes in all detected objects and generates new control signal
	}