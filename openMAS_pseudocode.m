
% Start File, ex: DSE_2_agents.m
% Set agent and obstacle type and count

% Run the simulation
environment/OMAS_initialize()
{
	% Set up simulation data and define agents and their public values
    % If you want to add any variables to META, add them to the struct in
    % this function
	ConfigureSimulationParameters
    
	% Simulate the system
	environment/OMAS_process()
	{
		% Simulate time steps
		ProcessGlobalTimeSeries()
		{
			for each time step:
			{
                % Step 0: Check agent idle status
                % Step 1: Update Timestep
                % Step 2: Record global state/trajectories
                % Step 3: Update the global states
                %   If you need to add noise to or access the true states
                %   do it here.
                environemnt/OMAS_updateGlobalStates()
                % Step 4: Update the visual representations
                %   IDK what this is
                % Step 5: Update communication matrix
                %   Call update obs and comm model functions
                % Step 6: Update simulation data
                %   Update collisions and relative positions
                % Step 7: Compute Information Filter Estimation
                %   Call a function to pass measurements to and run the
                %   information filter
                UpdateInformationFilter()
                {
                    % Compute observations from each agent to each other
                    % and pass those into the information filter, located
                    % in objects/agent.m, and inherited by every agent
                    referenceObject.InformationFilter()
                    {
                        % Compute one step of the information filter, and
                        % store priors and observations for consensus
                    }
                }
                % Step 8: Compute Consensus
                %   Break the agents into groups and re-sort agent states
                %   so all have the same ids in the same order.
                %   Run consensus
                consensus()
                % Step 9: Record Estimated State
                %   Record each agent's state estimates in the DATA struct
                %   so that it can be plotted later
                % Step 10: Update agent estimate from consensus data
                %   Give each agent the output of the consensus algorithm
                %   so it can update its controller. 
                % Step 11: Increment step counter
			}
		}
	}
    
    % After the simulation ends, generate plots
    OMAS_analysis()
}

Agents (In general)
Methods:
	main
	{
		Takes in all detected objects and generates new control signal
	}