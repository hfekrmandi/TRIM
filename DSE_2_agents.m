%% DEBUG SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
dbstop if error
rng(1);

% ADD THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
%addpath('toolboxes');
addpath('scenarios'); 

% CODE DEBUG COMMANDS %
% profile on
% profile viewer

fprintf('[SETUP]\tInitialising example script.\n');

%% INITIALISE ANY TOOLBOXES
% IntLab();   % Load Intlab
% OMAS_objectDiagnostics()

%% SIMULATION PARAMETERS
[~, userdir]   = system('echo %USERPROFILE%');  % Get desktop path
sim_outputPath = strcat(userdir,'\desktop\openmas-data');
sim_vebosity   = 1;
sim_warningDistance = 2;
sim_maxDuration = 30;
sim_timeStep    = 0.5;                        % Nominal (0.25s)
sim_idleTimeOut = 5*sim_timeStep;

sim_publishFigures = false;
% sim_publishFigures = true;
% sim_figureSet = {'all'};
% sim_figureSet = {'events','plan','inputs','isometric','gif'};
% sim_figureSet = {'plan','avoidance','inputs','isometric','gif'};
sim_figureSet = {'plan','estimates_rel', 'estimate_errors'};

%% SCENARIO PARAMETERS

%% Set obs and comm radii
% % Dead reckoning
% sim_obsRadius       = 0;
% sim_commRadius      = 0;

% % Information(/Kalman) Filter
% sim_obsRadius       = inf;
% sim_commRadius      = 0;

% Hybrid Consensus
sim_obsRadius       = inf;
sim_commRadius      = 0;

% Other parameters
sim_agentNumber     = 3;
sim_agentRadius     = 0.5;
sim_agentOrbit      = 3;
sim_agentVelocity   = 5;
sim_adjacencyMatrix = double(~eye(sim_agentNumber));
sim_waypointOrbit   = 2;
sim_waypointRadius  = 0.5;
sim_offsetAngle     = pi/4;
sim_obstacleNumber  = 0;
sim_obstacleOrbit   = 10;
sim_noiseSigma      = 0.2;
sim_plotScenario    = true;

%% INITIALISE AGENTS
fprintf('[SETUP]\tAssigning agent definitions:\n');
for index = 1:sim_agentNumber
% Listed here are the agent types. Uncomment one to use that agent in the
% simulation. 
    
% BASIC CLASSES
%     agentIndex{index} = objectDefinition('radius',sim_agentRadius);     
%     agentIndex{index} = agent('radius',sim_agentRadius);
%     agentIndex{index} = agent_test('radius',sim_agentRadius);
     agentIndex{index} = agent_2D();  
%     agentIndex{index} = agent_2D_test('radius',sim_agentRadius);

%     agentIndex{index} = agent_example('radius',sim_agentRadius);

% QUADCOPTER DYNAMICS
%     agentIndex{index} = quadcopter_legacy();
%     agentIndex{index} = quadcopter();
%     agentIndex{index} = quadcopter_formation('adjacencyMatrix',sim_adjacencyMatrix);
    
% ARdrone DYNAMICS
%     agentIndex{index} = ARdrone_prev();
%     agentIndex{index} = ARdrone('radius',sim_agentRadius);
%     agentIndex{index} = ARdrone_LQR();
%     agentIndex{index} = ARdrone_MPC();


%     agentIndex{index} = fixedWing();
%     agentIndex{index} = globalHawk();
%     agentIndex{index} = boeing737();
%     agentIndex{index} = A10();

%     agentIndex{index} = planetoid();
%     agentIndex{index} = earth();
%     agentIndex{index} = moon();
%     agentIndex{index} = ISS();

% FORMATION CONTROL 
%     agentIndex{index} = agent_formation('adjacencyMatrix',sim_adjacencyMatrix);
%     agentIndex{index} = agent_formation_boids();
%     agentIndex{index} = agent_formation_VO();
%     agentIndex{index} = agent_formation_RVO();
%     agentIndex{index} = agent_formation_HRVO();
%     agentIndex{index} = agent_2D_formation_VO('adjacencyMatrix',sim_adjacencyMatrix);
%     agentIndex{index} = agent_2D_formation_RVO();
%     agentIndex{index} = agent_2D_formation_HRVO();
%     agentIndex{index} = agent_2D_formation_RVO2();

% VECTOR SHARING;
%     agentIndex{index} = agent_vectorSharing('radius',sim_agentRadius,'detectionRadius',25);
%     agentIndex{index} = agent_2D_vectorSharing('radius',sim_agentRadius);

% INTERVAL AVOIDANCE
%     agentIndex{index} = agent_interval();
%     agentIndex{index} = agent_IA('radius',sim_agentRadius);
%     agentIndex{index} = agent_2D_IA('radius',sim_agentRadius);

% VELOCITY OBSTACLE METHODS
%     agentIndex{index} = agent_VO('radius',sim_agentRadius);
%     agentIndex{index} = agent_RVO('radius',sim_agentRadius);
%     agentIndex{index} = agent_HRVO('radius',sim_agentRadius);
%     agentIndex{index} = agent_2D_VO('radius',sim_agentRadius,'detectionRadius',25);
%     agentIndex{index} = agent_2D_unicycle('radius',sim_agentRadius,'detectionRadius',25);
%     agentIndex{index} = agent_2D_RVO('radius',sim_agentRadius);
%     agentIndex{index} = agent_2D_HRVO('radius',sim_agentRadius);
%     agentIndex{index} = agent_2D_RVO2('radius',sim_agentRadius); 
%     agentIndex{index} = agent_2D_VO_withComplex('radius',sim_agentRadius);
    
% OBSTACLES
%     agentIndex{index} = obstacle();
%     agentIndex{index} = obstacle_cuboid();
%     agentIndex{index} = obstacle_spheroid();
end

for index = 1:sim_obstacleNumber
	% OBSTACLES
%     obstacleIndex{index} = obstacle();
%     obstacleIndex{index} = obstacle_cuboid();
     %obstacleIndex{index} = obstacle_spheroid();
end

%% PLACE AGENT OBJECTS IN PRE-DEFINED SCENARIO
% FORMATION CONTROL TESTS
% [ objectIndex ] = GetScenario_corridor('agents',agentIndex,'adjacencyMatrix',sim_adjacencyMatrix,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_formation_split('agents',agentIndex,'agentSpacing',4,'adjacencyMatrix',sim_adjacencyMatrix,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = GetScenario_formation_fourObstacles('agents',agentIndex,'obstacles',sim_obstacleNumber,'obstacleRadius',1,'obstacle_orbit',sim_obstacleOrbit,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
 [ objectIndex ] = Scenario_3_agents_5_waypoints('agents',agentIndex,'obstacles',sim_obstacleNumber,'obstacleRadius',2,'obstacle_orbit',sim_obstacleOrbit,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma,'obsRadius',sim_obsRadius,'commRadius',sim_commRadius);

% % OBSTACLE TESTS
% [ objectIndex ] = GetScenario_fourCuboidObstacles('agents',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_obstacleTrack('agents',agentIndex,'obstacles',obstacleIndex,'plot',sim_plotScenario);

% % AGENT TESTS
% [ objectIndex ] = GetScenario_twoLines('agents',agentIndex,'agentVelocity',sim_agentVelocity,'padding',5,'agentSeparation',sim_agentOrbit,'waypointSeparation',sim_waypointOrbit,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'waypointRadius',sim_waypointRadius,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = GetScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = GetScenario_concentricSphere('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_concentricAngle('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'angle',sim_offsetAngle,'plot',sim_plotScenario);

% RANDOM TESTS
% [ objectIndex ] = GetScenario_random('objects',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_randomNormal('objects',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_randomUniform('objects',agentIndex,'is3D',false,'plot',sim_plotScenario);

% WAYPOINT TESTS
% [ objectIndex ] = GetScenario_waypointCurve('agents',agentIndex,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_waypoint_90Degrees('agents',agentIndex,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);

% % RL TESTS
% [ objectIndex ] = GetScenario_earthOrbit('plot',sim_plotScenario);

%% %%%%%% INITIALISE THE SIMULATION WITH THE OBJECT INDEX %%%%%%%%%%%%%%%%%
[DATA,META] = OMAS_initialise('objects',objectIndex,...
                             'duration',sim_maxDuration,...
                                   'dt',sim_timeStep,...
                          'idleTimeOut',sim_idleTimeOut,...
                              'figures',sim_figureSet,...
                      'warningDistance',sim_warningDistance,...
                            'verbosity',sim_vebosity,...
                           'outputPath',sim_outputPath,...
                          'publishMode',sim_publishFigures);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except DATA META
load(strcat(META.outputPath,'META.mat'));
load(strcat(META.outputPath,'EVENTS.mat'));
load(strcat(META.outputPath,'OBJECTS.mat'));

%% EXTRA GRAPHS (DEBUG)
% objectDATA = DATA.objectIndex{1,1}.DATA;
% figure()
% grid on; 
% plot(META.TIME.timeVector,objectDATA.inputs)
% xlabel('time (s)');
% legend('Location','southeast');
