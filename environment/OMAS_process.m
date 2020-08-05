%% THE OPENMAS CYCLE PROCESSOR (OMAS_process.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the progressive cycles of the simulation, for the
% environment, obstacles, agents and objects withing the simulation space.

% Author: James A. Douthwaite 06/10/2016

% //////////////////////// MAIN WRAPPER SCRIPT ////////////////////////////
function [DATA,SIM,EVENTS,objectIndex]  = OMAS_process(META,objectIndex)
% INPUTS:
% objectIndex - The cell array of object classes
% OUTPUTS:
% DATA        - The output DATA structure with timeseries data
% SIM         - The terminal META structure copy

%% PRE-SIMULATION DECLARATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'META'; % Declare the simuation active phase
EVENTS = [];         % Reset the EVENT log container

% PREPARE THE OUTPUT DATA CONTAINER
[DATA] = GetOutputStructure(META);

%%%%%%%%%%%%%%%%%%%%%% BEGIN TIME-STEP INTERATIONS %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\n[%s]\tLAUNCHING SIMULATION...\n[%s]\n',META.phase,META.phase,META.phase);
% EXECUTE OMAS MAIN CYCLE
[META,objectIndex,DATA,EVENTS] = ProcessGlobalTimeSeries(META,objectIndex,DATA,EVENTS);
% DATA.computationTime = toc; % Get total elapsed simulation time
% fprintf('\n[%s]\tSIMULATION COMPLETE (time elapsed: %ss)\n',META.phase,num2str(DATA.computationTime));
fprintf('[%s]\n[%s]\t...SIMULATION COMPLETE\n[%s]\n',META.phase,META.phase,META.phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINALISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'OUTPUT';                                   % Change to OUTPUT phase
% clearvars -except META EVENTS DATA objectIndex         % Clear loose simulation variables
fprintf('[%s]\tDumping initial simulation variables to output directory\n[%s]\n',META.phase,META.phase);
% sendOutputToFiles(META,EVENTS,DATA,objectIndex);
fprintf('[%s]\tClosing simulation.\n',META.phase);
SIM = META;
end

% //////////////////////// MAIN CYCLE FUNCTIONS ///////////////////////////
% PROCESS TIME VECTOR
function [META,objectIndex,DATA,EVENTS] = ProcessGlobalTimeSeries(META,objectIndex,DATA,EVENTS)
% This function computes the simulation loop cycle across the given
% descrete META.timevector:
% INPUTS:
% META - The initialised META data-structure
% objectIndex - The cell vector of objects (and/or agents)
% DATA        - The initialised output DATA structure
% EVENTS      - The initialised EVENTS structure list

% PROCESS TIME CYCLE
step = META.TIME.currentStep; 
exitFlag = 0; 
while step <= META.TIME.numSteps
    
    % 0. /////////////////// CHECK AGENT IDLE STATUS //////////////////////
    % If there are agents in the simulation, then break if all are idle,
    % else allow run to full duration
    idleCondition = sum([META.OBJECTS([META.OBJECTS.type] == OMAS_objectType.agent).idleStatus]) == META.totalAgents && META.totalAgents > 0;
    if idleCondition == 1 && ~exitFlag
        exitFlag = 1; 
        referenceTime = META.TIME.currentTime;
        fprintf('[%s]\t... All agents are idle, aborting in t=%0.2fs.\n',META.phase,META.TIME.idleTimeOut);
    elseif exitFlag == 1
        timeOut = META.TIME.idleTimeOut - (META.TIME.currentTime - referenceTime);
        if timeOut <= 0
           fprintf('[%s]\t... Aborting.\n',META.phase); 
           break 
        end        
        fprintf('[%s]\t... All agents are idle, aborting in t=%0.2fs.\n',META.phase,timeOut);
    end 
    % /////////////////////////////////////////////////////////////////////
    
    % 1. ////////////// UPDATE TIMESTEP PARAMETERS (@t=k) /////////////////
    META.TIME.currentTime = META.TIME.timeVector(step);                    % The current sim-time
    META.TIME.currentStep = step;                                          % The current sim-step
    if META.verbosity >= 1
        fprintf('[%s]\tStep: %.0f\tTime: %.2fs\n',META.phase,META.TIME.currentStep,META.TIME.currentTime);
    end  
    % /////////////////////////////////////////////////////////////////////
    
    % 2. ///////////////// RECORD GLOBAL STATE (@t=k) /////////////////////
    for ID1 = 1:META.totalObjects
        % Collect the META.OBJECTS.state data (as the simulations understanding of the global states) 
        % and save to the output DATA.globalTrajectories set, must be done synchronously).
        DATA.globalTrajectories(DATA.stateIndex(1,ID1):DATA.stateIndex(2,ID1),META.TIME.currentStep) = META.OBJECTS(ID1).globalState;
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 3. ////////////// UPDATE THE GLOBAL STATES (@t=k) ////////////////////
    for ID1 = 1:META.totalObjects                                          % For each META object element
        % Update META.OBJECT the structure
%         META.OBJECTS(ID1) = OMAS_updateGlobalStates_mex(...
%             META,...
%             objectIndex{ID1}.objectID,...
%             objectIndex{ID1}.GetGLOBAL('velocity'),...
%             objectIndex{ID1}.GetGLOBAL('quaternion'),...
%             objectIndex{ID1}.GetGLOBAL('idleStatus'));   
        META.OBJECTS(ID1) = OMAS_updateGlobalStates(...
            META,...
            objectIndex{ID1}.objectID,...
            objectIndex{ID1}.GetGLOBAL('velocity'),...
            objectIndex{ID1}.GetGLOBAL('quaternion'),...
            objectIndex{ID1}.GetGLOBAL('X'),...
            objectIndex{ID1}.GetGLOBAL('idleStatus')); 
    end
    % /////////////////////////////////////////////////////////////////////
    
%     % 4. ////////// UPDATE THE VISUAL REPRESENTATIONS (@t=k) //////////////
%     for ID1 = 1:META.totalObjects
%         if ~isstruct(META.OBJECTS(ID1).geometry)
%             continue
%         end
%         % MODIFY THE ASSOCIATED GEOMETRY
%         META.OBJECTS(ID1) = updateVisuals(META.OBJECTS(ID1));
%     end
    
    %  5. ////////////// UPDATE Communication matrix (@t=k) ///////////////
%     if step < 20 || step > 40
%         comm_list = {{2, 3}; {1, 3}; {1, 2}};
%     else
%         comm_list = {{3}; {3}; {1, 2}};
%     end
%     comm_list = {{}; {}; {}};
    objectIndex = apply_obs_model_distance(META, objectIndex);
    objectIndex = apply_comm_model_distance(META, objectIndex);

    % 6. //////// UPDATE SIMULATION/ENVIRONMENTAL META DATA (@t=k) ////////
    [META,metaEVENTS] = UpdateSimulationMeta(META,objectIndex);            % Update META snapshot with equivilant objectIndex.state
    % LOG THE META EVENTS
    if ~isempty(metaEVENTS)                                                % If META events occurred this timestep
        EVENTS = vertcat(EVENTS,metaEVENTS);                               % Append to global EVENTS
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 7. ////////// COMPUTE INFORMATION FILTER ESTIMATION (@t=k) //////////
    if META.threadPool ~= 0
        objectSnapshot = objectIndex;                                      % Make a temporary record of the object set
        parfor (ID1 = 1:META.totalObjects)
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [detection{ID1},objectIndex{ID1},objectEVENTS] = UpdateInformationFilter(META,objectSnapshot,objectIndex{ID1});
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)
                EVENTS = vertcat(EVENTS,objectEVENTS);
            end
        end
    else
        for ID1 = 1:META.totalObjects
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [detection{ID1},objectIndex{ID1},objectEVENTS] = UpdateInformationFilter(META,objectIndex,objectIndex{ID1}); % Update objectIndex snapshot with new META data
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)                                      % If objectEVENTS occur in this timestep
                EVENTS = vertcat(EVENTS,objectEVENTS);                     % Append to global EVENTS
            end
        end
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 8. ///////////////// COMPUTE CONSENSUS STEPS (@t=k) /////////////////
    agent_data = get_sorted_agent_states(META, objectIndex);
    agent_groups = break_agents_into_groups(META, agent_data);
    consensus(agent_groups);
    % /////////////////////////////////////////////////////////////////////
    
    % 9. ///////////////// RECORD ESTIMATED STATE (@t=k) /////////////////////
    for index = 1:numel(agent_data)
        % Collect the META.OBJECTS.state data (as the simulations understanding of the global states) 
        % and save to the output DATA.globalTrajectories set, must be done synchronously).
        agent = agent_data{index};
        X_agent = agent.memory_x;
        P_agent = agent.memory_P;
        IDs = agent.memory_id_list;
        Obs = agent.memory_id_obs;
        Comm = agent.memory_id_comm;
        DATA.ids(index, 1:numel(IDs), META.TIME.currentStep) = IDs';
        DATA.Observations(index, 1:numel(Obs), META.TIME.currentStep) = Obs';
        DATA.Connections(index, 1:numel(Comm), META.TIME.currentStep) = Comm';
        for ID = IDs
            P = agent.cov_from_id(P_agent, IDs, agent.objectID);
            DATA.estimate_covariance(index, ID, :, :, META.TIME.currentStep) = P;
            
            X_zero = agent.state_from_id(X_agent, IDs, agent.objectID);
            X_estimate = agent.state_from_id(X_agent, IDs, ID);
            X_relative = X_estimate - X_zero;
            DATA.estimates_rel(index, ID, :, META.TIME.currentStep) = X_relative;

            X_gt_agent = META.OBJECTS(agent.objectID).X;
            X_gt_object = META.OBJECTS(ID).X;
            X_gt_relative = X_gt_object - X_gt_agent;
            DATA.estimate_errors(index, ID, :, META.TIME.currentStep) = X_relative - X_gt_relative;
        end
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 10. //////// UPDATE AGENT ESTIMATE FROM CONSENSUS DATA (@t=k) ////////
    if META.threadPool ~= 0
        objectSnapshot = objectIndex;                                      % Make a temporary record of the object set
        parfor (ID1 = 1:META.totalObjects)
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [objectIndex{ID1},objectEVENTS] = UpdateObject(META,objectSnapshot,objectIndex{ID1});            
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)
                EVENTS = vertcat(EVENTS,objectEVENTS);
            end
        end
    else
        for ID1 = 1:META.totalObjects
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [objectIndex{ID1},objectEVENTS] = UpdateObject(META,objectIndex,objectIndex{ID1}); % Update objectIndex snapshot with new META data
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)                                      % If objectEVENTS occur in this timestep
                EVENTS = vertcat(EVENTS,objectEVENTS);                     % Append to global EVENTS
            end
        end
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 11. /// THE 'OBJECT.VIRTUAL' PROPERTIES IS NOW UPDATED FOR (t=k+1) ///
    step = step + 1;
end
% CREATE TERMINAL VALUES, FOR CLARITY
META.TIME.endStep = META.TIME.currentStep;
META.TIME.endTime = META.TIME.currentTime; 
end
% UPDATE THE OBJECT META PROPERTIES
function [SIM,metaEVENTS]	= UpdateSimulationMeta(SIM,objectIndex)
% In this function we wish to update the global separations and event
% conditions based on the new 'globalState' properties.

% INPUTS:
% SIM         - A local copy of the global META structure
% objectIndex - The cell array of object classes.
% OUTPUTS:
% SIM         - The updated META structure
% metaEVENTS  - A vector of new event objects.

% We want to move through the object set an update both the object being
% updated, but also the object it is being updated against.

% //// ASSESS THE COLLISIONS AND RELATIVE POSITIONS (SYMMETRICAL CHECK) ///
collisionLogicals = zeros(SIM.totalObjects);
warningLogicals   = zeros(SIM.totalObjects);
for entityA = 1:SIM.totalObjects                                           % Object A's position in the META.OBJECTS
    
    % NOTE:
    % - We only need to consider the objects whos separations have yet to be
    %   evaluated. If we update 1&2 and 2&1 simultaneously, then we must only
    %   examine the unique ID permutations to assess their conditions.
    
    for entityB = (entityA+1):(SIM.totalObjects)                           % Object B's position in the META.OBJECTS
        % THE OBJECTS ASSOCIATED objectIndex
        object_A = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityA).objectID};
        object_B = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID};  % Determine the associated geometries
        
        % UPDATE RELATIVE POSITIONS  (centroid separation based)
        SIM.OBJECTS(entityA).relativePositions(entityB,:) = (SIM.OBJECTS(entityB).globalState(1:3,1) - SIM.OBJECTS(entityA).globalState(1:3,1))';   % AB vector
        SIM.OBJECTS(entityB).relativePositions(entityA,:) = (SIM.OBJECTS(entityA).globalState(1:3,1) - SIM.OBJECTS(entityB).globalState(1:3,1))';   % BA vector
        
        % GET THE WARNING CONDITIONS (centroid separation based)
        separationDistance = norm(SIM.OBJECTS(entityA).relativePositions(entityB,:)) - (SIM.OBJECTS(entityA).radius + SIM.OBJECTS(entityB).radius);
        warningCondition = SIM.warningDistance >= separationDistance;
        warningLogicals(entityA,entityB) = warningCondition;
        warningLogicals(entityB,entityA) = warningCondition;               % Assign the warning condition
        
        % EVALUATE COLLISIONS BETWEEN THE OBJECTS
        collisionCondition = OMAS_collisionDetection(...
            SIM.OBJECTS(entityA),object_A.GEOMETRY,...
            SIM.OBJECTS(entityB),object_B.GEOMETRY,...
            SIM.conditionTolerance);
        collisionLogicals(entityA,entityB) = collisionCondition;
        collisionLogicals(entityB,entityA) = collisionCondition;           % Assign the collision condition     
        
        % ERROR CHECKING (Should never happen)
        if separationDistance > 0 && collisionCondition
            error('[ERROR] A collision occurred at distance %f without violating the minimum separation of %f',...
                    norm(SIM.OBJECTS(entityA).relativePositions(entityB,:)),(SIM.OBJECTS(entityA).radius + SIM.OBJECTS(entityB).radius));
        end
    end
end

% //////////// ASSESS DETECTION CONDITIONS (ASYMMETRICAL CHECK) ///////////
% DEFAULT DETECTION CONDITION
detectionLogicals = zeros(SIM.totalObjects);
for entityA = 1:SIM.totalObjects
    % Check objects have the capacity to observe other objects
    if SIM.OBJECTS(entityA).type ~= OMAS_objectType.agent
        continue
    end
    
    % Check all agent observations
    for entityB = 1:SIM.totalObjects                                       % Object B's position in the META.OBJECTS
        % SKIP SELF-REFERENCE    
        if SIM.OBJECTS(entityA).objectID == SIM.OBJECTS(entityB).objectID                                           
            continue                                                       % Only agents can generate notifications
        end 
        % NOTE:
        % - We also need to determine which agents are able to observe other
        %   agents in the field. However, this is NOT A SYMMETRIC operation
        %   as A may be able to detect B without B detecting A. 
        
        % Detection check #1 - Authorisation
        isWaypoint = SIM.OBJECTS(entityB).type == OMAS_objectType.waypoint;
        if isWaypoint
            isAuthorised = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID}.IDAssociationCheck(SIM.OBJECTS(entityA).objectID);
            if ~isAuthorised
                continue
            end
        end
        
        % Detection check #2 - Spherical 
        isDetected = OMAS_geometry.intersect_spheres(...
            SIM.OBJECTS(entityA).globalState(1:3),SIM.OBJECTS(entityA).detectionRadius + 0.5*SIM.conditionTolerance,...
            SIM.OBJECTS(entityB).globalState(1:3),SIM.OBJECTS(entityB).radius + 0.5*SIM.conditionTolerance);
        
        % Object is not detected
        if ~isDetected
            continue
        end
        
        % This check is sufficient for way-points
        if isWaypoint
            detectionLogicals(entityA,entityB) = 1;
            continue
        end
        
        % ///////////////// OBJECT MAY BE OBSERVED ////////////////////////
        % NOTE:
        % - If the distance(less the radius) is less than detection radius
        %   then its worth checking for more a more complex representation.
        % - We want to check if any vertices/edges can be observed. If any 
        %   edges are visible, then they are to be sent to the first object.
        
        % Get the geometry of the second object (GEOMETRY OF B)
        geometryB = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID}.GEOMETRY;        
         
        % Check if a better check is necessary
        hasGeometry = isstruct(geometryB) && size(geometryB.vertices,1) > 0;    
        if isDetected && ~hasGeometry 
            detectionLogicals(entityA,entityB) = 1;                        % Detection occurred
            break
        end
        isDetected = 0; % Reset, potential false positive when checking against polygon
         
        % Assess the polygon
        relativeR = SIM.OBJECTS(entityA).R'*SIM.OBJECTS(entityB).R;        % Rotate second geometry orientated relative to
        % Design the vertex set
        geometryB.vertices = geometryB.vertices*relativeR + SIM.OBJECTS(entityB).globalState(1:3)';  
        
        % Assess each face and vertex for detection
        for face = 1:size(geometryB.vertices,1)
            % GET THE MEMBER IDs OF THE FACE
            faceMembers  = geometryB.faces(face,:);
            faceVertices = geometryB.vertices(faceMembers,:);
            
            % CHECK THE GEOMETRY AGAINST THE SPHERICAL CONSTRAINT (with mex)
%             isDetected = CheckSphereTriangleIntersection_mex(...
%                 SIM.OBJECTS(entityA).globalState(1:3),...
%                 SIM.OBJECTS(entityA).detectionRadius,...
%                 faceVertices(1,:)',faceVertices(2,:)',faceVertices(3,:)');
            % CHECK THE GEOMETRY AGAINST THE SPHERICAL CONSTRAINT
            isDetected = CheckSphereTriangleIntersection(...
                SIM.OBJECTS(entityA).globalState(1:3),...
                SIM.OBJECTS(entityA).detectionRadius,...
                faceVertices(1,:)',faceVertices(2,:)',faceVertices(3,:)');
            
            % IF THE FACE IS DETECTED
            if isDetected
                detectionLogicals(entityA,entityB) = 1;
                break
            end
        end
        
        % NOTE:
        % - If any faces are within the constraint, then the detection
        %   logical is set true and then exited.
        % - Otherwise, if the none of the contraints are met, then it
        %   remains undetected.
    end
end

% NOTES:
% - The fields geometry based fields of SIM.OBJECTS are updated to this
%   time-step. The event conditions can now be re-evaluted.
% - CHECKS MUST BE ASYMMETRICAL - A may detect B without B detecting A

metaEVENTS = [];                                                           % Reset meta Events container
for entityA = 1:SIM.totalObjects
    % CHECK OBJECTS HAVE THE CAPACITY TO OBSERVE OTHER OBJECTS
    if SIM.OBJECTS(entityA).type ~= OMAS_objectType.agent
        continue
    end
    
    for entityB = 1:SIM.totalObjects
        % OMIT SELF-CHECK
        if SIM.OBJECTS(entityA).objectID == SIM.OBJECTS(entityB).objectID                                         
            continue                                                       % Only agents can generate notifications
        end 
                
        % //////////// ASSESS EVENT CONDITIONS FOR THE AGENTS /////////////
        % NOTES:
        % - The agent set is only capable of generating events associated
        %   with the detections, warnings and collisions.
        
        % ///////////////// ASSESS DETECTIONS CONDITIONS //////////////////
        % DETECTION EVENT LOGIC
        ABdetectionConstraint       =  detectionLogicals(entityA,entityB);
        novelDetectionConstraint    =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.detection) == 0;
        detectionEventCondition     =  ABdetectionConstraint &&  novelDetectionConstraint;   % Is within a range and a detection event has not been issued.
        detectionEventNullCondition = ~ABdetectionConstraint && ~novelDetectionConstraint;   % Is outside a range and a detection condition is still toggled.
        % EVALUATE OCCURANCE
        if detectionEventCondition 
            % UPDATE DETECTION STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.detection) = 1;              % Ammend status of the META object
            % GENERATE DETECTION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.detection);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        elseif detectionEventNullCondition
            % AMEND THE META OBJECT
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.detection) = 0;              % Ammend status of the META object
            % GENERATE THE DETECTION NULLIFICATION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_detection);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        end
                
        % ///////////////// ASSESS WAYPOINT CONDITIONS ////////////////////
        if SIM.OBJECTS(entityB).type == OMAS_objectType.waypoint
            % Check A is allowed to get B
            waypointObj = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID};
            hasAssociation = waypointObj.IDAssociationCheck(SIM.OBJECTS(entityA).objectID);
            % WAYPOINT EVENT LOGIC
            ABwaypointCondition        =  collisionLogicals(entityA,entityB) && hasAssociation;
            novelWaypointGetConstraint =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.waypoint) == 0;
            waypointGetEventCondition  =  ABwaypointCondition &&  novelWaypointGetConstraint;    % Is within a range and a collision event has not been issued.
            waypointEventNullCondition = ~ABwaypointCondition && ~novelWaypointGetConstraint;    % Is outside a range and a waypoint-get condition is still toggled.
            % EVALUATE OCCURANCE
            if waypointGetEventCondition
                % AMEND THE META OBJECT
                SIM.OBJECTS(entityA).objectStatus(entityB,eventType.waypoint) = 1;
                % GENERATE THE WAYPOINT ACHIEVED EVENT
                [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.waypoint);
                metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
            elseif waypointEventNullCondition
                % AMEND THE META OBJECT
                SIM.OBJECTS(entityA).objectStatus(entityB,eventType.waypoint) = 0;
                % GENERATE THE WAYPOINT NULLIFICATION EVENT
                %[EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_waypoint);
                %metaEVENTS = vertcat(metaEVENTS,EVENT);                                        % Add new META EVENTS to global structure
            end 
            continue    % Skip collision check
        end

        % ////////// CONFIRM THE OBJECT CAN BE COLLIDED WITH //////////////      
        if SIM.OBJECTS(entityB).hitBox == OMAS_hitBoxType.none
            % Check hit-box definition:
            % -If no hitbox is assigned, no need to evaluate warnings or
            % collisions.
            continue              
        end
        
        % /////////////////// ASSESS WARNING CONDITIONS ///////////////////
        % WARNING EVENT LOGIC
        ABwarningConstraint       =  warningLogicals(entityA,entityB);
        novelWarningConstraint    =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.warning) == 0;
        warningEventCondition     =  ABwarningConstraint &&  novelWarningConstraint;         % Is within a range and a warning event has not been issued.
        warningEventNullCondition = ~ABwarningConstraint && ~novelWarningConstraint;         % Is outside a range and a warning condition is still toggled.       
        % EVALUATE OCCURANCE
        if warningEventCondition 
            % UPDATE WARNING STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.warning) = 1;
            % GENERATE THE WARNING EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.warning);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        elseif warningEventNullCondition
            % UPDATE WARNING-NULL STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.warning) = 0;
            % GENERATE THE WARNING NULLIFICATION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_warning);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        end
        
        % ////////////////// ASSESS COLLISION CONDITIONS //////////////////
        % COLLISION EVENT LOGIC
        ABcollisionConstraint       =  collisionLogicals(entityA,entityB);
        novelCollisionConstraint    =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.collision) == 0;   
        collisionEventCondition     =  ABcollisionConstraint &&  novelCollisionConstraint;   % Is within a range and a collision event has not been issued.
        collisionEventNullCondition = ~ABcollisionConstraint && ~novelCollisionConstraint;   % Is outside a range and a collision condition is still toggled.     
        % EVALUATE OCCURANCE
        if collisionEventCondition
            % UPDATE COLLISIONL STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.collision) = 1;
            % GENERATE THE COLLISION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.collision);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        elseif collisionEventNullCondition
        	% AMEND THE META OBJECT
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.collision) = 0;
            % GENERATE THE COLLISION NULLIFICATION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_collision);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        end     
    end
end

end
% UPDATE THE OBJECT PROPERTIES
function [observationPacket, referenceObject,objectEVENTS] = UpdateInformationFilter(SIM,objectIndex,referenceObject)
% This function updates a referenceObject class against the rest of the
% objectIndex, independantly of the SIM.OBJECTS META data.
% INPUTS:
% SIM             - Current META structure
% objectIndex     - The current object list
% referenceObject - The agent object being updated

% OUTPUTS:
% referenceObject - The updated agent class

objectEVENTS = []; % Container for object based events

% GET THE OBJECTS EQUIVALENT SIM OBJECT
SIMfirstObject = SIM.OBJECTS(SIM.globalIDvector == referenceObject.objectID);

% If the observing object is "2D"
if referenceObject.Is3D()
    dimensionIndices = 1:3;
else
    dimensionIndices = 1:2;
end

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling ID:%d\t name: %s\n',SIMfirstObject.objectID,referenceObject.name);

end

% Default timing parameters
ENV = SIM.TIME;
ENV.outputPath = SIM.outputPath;

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling ID:%d\t name: %s\n',SIMfirstObject.objectID,referenceObject.name);
end

observationPacket = [];
% SWITCH BEHAVIOUR BASED ON OBJECT TYPE
switch SIMfirstObject.type
    case OMAS_objectType.agent
        % AGENT - SIMULATION/ENVIROMENTAL FEEDBACK REQUIRED %%%%%%%%%%%%%%
        % Container for agent detection packet
        
        % //// PARSE THE ENVIRONMENTAL OBJECTS OBSERVABLE TO THE AGENT ////
        for ID2 = 1:SIM.totalObjects
            % GET THE SIM OBJECT OF THE EVALUATION OBJECT
            SIMsecondObject = SIM.OBJECTS(ID2);
            firstObject = objectIndex{SIMfirstObject.objectID};
            
            isSameObject = SIMfirstObject.objectID == SIMsecondObject.objectID;
            
            obs_list = firstObject.memory_id_obs;
            index = find(obs_list == ID2);
            if ~isempty(index)
                isObserved = true;
            else
                isObserved = false;
            end
            
%             % Skip condition - Get the current status of this object
%             isDetected = 1 == SIMfirstObject.objectStatus(ID2,eventType.detection);
%             isSameObject = SIMfirstObject.objectID == SIMsecondObject.objectID; 
%             if ~isDetected || isSameObject
%                 continue
%             end
            
            if SIMsecondObject.type == OMAS_objectType.waypoint
                warning('off')
                association = objectIndex{ID2}.GetAgentAssociation(referenceObject);
                if numel(association) < 1
                    continue
                end
                warning('on')
            end

            if (~isObserved || isSameObject) && SIMsecondObject.type == 1
                continue
            end
            
            % GET THE EQUIVALENT OBJECT TO EVALUATE AGAINST THE REFERENCE
            secondObject = objectIndex{SIM.globalIDvector == SIMsecondObject.objectID};
            
            % ///////////////// ELSE AGENT IS DETECTED ////////////////////
            % ////// (GENERATE A LOCALISED PROXIMITY DESCIPTION) //////////
            % The second agent is within global seperation conditions. A
            % communication packet is therefore assembled containing the
            % data on the second object. rotated into the agent-local
            % coordinate frame.

            % GET THE RELATIVE POSITIONS OF THE AGENT AND OBJECT IN THE 
            % GLOBAL COORDINATE SYSTEM (x_rel = x_B - x_A)
            relativePosition = SIMfirstObject.relativePositions(ID2,:)';
            relativeVelocity = SIMsecondObject.globalState(4:6) - SIMfirstObject.globalState(4:6);
            % ROTATE THE GLOBAL STATE OF THE OBJECT INTO THE AGENT FRAME
            %observedPosition = SIMfirstObject.R*relativePosition;         % Rotate the from the global into the body frame of the simReference
            %observedVelocity = SIMfirstObject.R*relativeVelocity;
            observedPosition = relativePosition(1:2);                           % Rotate the from the global into the body frame of the simReference
            observedVelocity = relativeVelocity;
            observedTheta = SIMsecondObject.eulZYX(1);
            % SPHERICAL REPRESENTATION
            observedRange    = norm(observedPosition);
            
            %noise_sigma = 0.2;
            %noise = noise_sigma*eye(6);
            % Add camera-specific noise
            % [x y z phi theta psi dx dy dz dphi dtheta dpsi]
            relative_X = SIMsecondObject.X(1:6) - SIMfirstObject.X(1:6);
            range = norm(relative_X(1:3));
            noise = svgs_R_from_range_SRT(range);
            observedz = relative_X + mvnrnd(zeros(1, 6), noise)';
            
            %observedElevation = asin(observedPosition(3)/observedRange);
            %observedHeading   = atan2(observedPosition(2),observedPosition(1));
            % PRESENT SIZE PROPERTIES
            %observedRadius    = SIMfirstObject.radius;     
            %observedAngularWidth = 2*asin(observedRadius/(observedRange + observedRadius));               
                        
            % OBJECT GEOMETRY PREPARATION
            % In order for the agent to observe the object correctly. The
            % geometry must be translated and rotated into the relative
            % frame of the detecting agent.
            % ASSUMPTION:
            % - The geometry is normalised to the body axes of the object.
            % TO DO:
            % - If the object is on the edge of vision.. only part of the
            %   object will be visable to the agent. a process must be in
            %   place to create a subset of the geometry.

            % OBJECT PRIORITY (IF WAYPOINT)
            observedPriority = NaN;
            if SIMsecondObject.type == OMAS_objectType.waypoint
                association = objectIndex{ID2}.GetAgentAssociation(referenceObject);
                observedPriority = association.priority;
            end
            
            % A DEBUG STUCTURE 
            DEBUG = struct('globalPosition',SIMsecondObject.globalState(dimensionIndices),...     % Added for simplicity
                           'globalVelocity',SIMsecondObject.globalState(dimensionIndices + 3),... % Added for simplicity
                           'priority',observedPriority);                                          % The objects global priority for that agent
            
            % ///////////////// PREPARE DETECTION PACKET //////////////////
            detectionObject = struct('objectID',SIMsecondObject.objectID,...            % The object ID (observed)
                                     'name',SIMsecondObject.name,...                    % The object name tag (observed)
                                     'type',SIMsecondObject.type,...                    % The objects sim-type enum
                                     'z',observedz,...                                  % The apparent range
                                     'colour',SIMsecondObject.colour,...                % Finally the simulation's colourID 
                                     'DEBUG',DEBUG);                                    % Pass additional parameters (not otherwise known)
            observationPacket = vertcat(observationPacket,detectionObject);             % Append object to packet to agent
        end

        % /////////////// COMPUTE THE LOCAL AGENT CYCLE ///////////////////
        % Given the detection object defined for this agent, compute the
        % agents cycle with this information.
        % SEND OBJECT OBSERVERATION PACKET TO AGENT
        referenceObject.InformationFilter(ENV.dt,observationPacket);
    case OMAS_objectType.obstacle
        % PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        referenceObject = referenceObject.main(ENV);
    case OMAS_objectType.waypoint
        % WAYPOINTS ARE CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%
        referenceObject = referenceObject.main(ENV);
    otherwise
        % CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%
        % DO NOT UPDATE
end
end

% /////////////////////// CONSENSUS OPERATIONS ////////////////////////////

% Grab the agent with a specific ID
% Done by looping through the list of agents and checking the ID
function [agent] = agent_with_id(agents, ID)
    for i = 1:numel(agents)
        if agents(i).objectID == ID
            agent = agents(i);
            return;
        end
    end
    agent = [];
end

% Ensures that every agent has the same state variables in the same order
function [agents] = get_sorted_agent_states(SIM,objectIndex)

    % Build combined list of ids
    id_list = [];
    for agent_index = 1:numel(SIM.OBJECTS)
        if SIM.OBJECTS(agent_index).type == OMAS_objectType.agent
            for i = 1:numel(SIM.OBJECTS)
                if objectIndex{agent_index}.objectID == SIM.OBJECTS(agent_index).objectID
                    agents{objectIndex{agent_index}.objectID} = objectIndex{agent_index};
                    id_list = [id_list, objectIndex{agent_index}.memory_id_list];
                end
            end
        end
    end
    
    % Ensure that the list is sorted, so it is the same on sequential runs
    id_list = sort(unique(id_list));
    dim_state = agents{1}.dim_state;
    dim_obs = agents{1}.dim_obs;
    n_agents = numel(id_list);
    
    % Ensure all agents' state variables match the master list
    for agent_index = 1:numel(agents)
        agent = agents{agent_index};
        
        % If the state variables don't match, add them in
        if ~isequal(agent.memory_id_list, id_list)
            Y = 0.01*eye(n_agents*dim_state);
            y = zeros(n_agents*dim_state, 1);
            I = zeros(n_agents*dim_state);
            i = zeros(n_agents*dim_state, 1);
            
            % Move the agents' values to the location specified in the master list
            for agent_index_1 = 1:numel(agent.memory_id_list)
                for agent_index_2 = 1:numel(agent.memory_id_list)
                    
                    group_index_1 = find(id_list == agent.memory_id_list(agent_index_1));
                    group_index_2 = find(id_list == agent.memory_id_list(agent_index_2));
                    
                    % Generate indices (to make the assignment setp shorter)
                    g_row_lo = dim_state*(group_index_1 - 1) + 1;
                    g_row_hi = dim_state*group_index_1;
                    g_col_lo = dim_state*(group_index_2 - 1) + 1;
                    g_col_hi = dim_state*group_index_2;
                    a_row_lo = dim_state*(agent_index_1 - 1) + 1;
                    a_row_hi = dim_state*agent_index_1;
                    a_col_lo = dim_state*(agent_index_2 - 1) + 1;
                    a_col_hi = dim_state*agent_index_2;
                    
                    Y(g_row_lo:g_row_hi,g_col_lo:g_col_hi) = agent.memory_Y(a_row_lo:a_row_hi,a_col_lo:a_col_hi);
                    I(g_row_lo:g_row_hi,g_col_lo:g_col_hi) = agent.memory_I(a_row_lo:a_row_hi,a_col_lo:a_col_hi);
                end
                
                y(g_row_lo:g_row_hi) = agent.memory_y(a_row_lo:a_row_hi);
                i(g_row_lo:g_row_hi) = agent.memory_i(a_row_lo:a_row_hi);
            end
            
            agent.memory_id_list = id_list;
            agent.memory_Y = Y;
            agent.memory_y = y;
            agent.memory_I = I;
            agent.memory_i = i;
        end
    end
end

% Update the momory_id_comm list in each agent to include each agent they
% can communicate with (based on a communication model).
function [agents] = apply_comm_model_obs(agents)

    % Apply communication model and create list of agents each agent can communicate with
    % Current communication model is the same as the observation model
    agents_arr = [];
    for agent = agents
        agents_arr = [agents_arr, agent{1}];
    end
    for agent = agents_arr
        agent.memory_id_comm = [];
        for i = 1:numel(agent.memory_id_obs)
            if ~isempty(agent_with_id(agents_arr, agent.memory_id_obs(i)))
                agent.memory_id_comm = [agent.memory_id_comm, agent.memory_id_obs(i)];
            end
        end
    end
end

% Update the momory_id_comm list in each agent to include each agent they
% can communicate with (based on a communication model).
function [agents] = apply_comm_model(SIM, agents, comm_list)

    % Apply communication model and create list of agents each agent can communicate with
    % Current communication model is the same as the observation model
    agents_arr = [];
    for id = 1:numel(agents)
        agent = agents{id};
        if SIM.OBJECTS(agent.objectID).type == 1
            agents_arr = [agents_arr, agent];
        end
    end
    for index = 1:numel(agents_arr)
        agent = agents_arr(index);
        agent_comm_list = cell2mat(comm_list{index});
        agent.memory_id_comm = agent_comm_list;
    end
end

% Update the momory_id_comm list in each agent to include each agent they
% can communicate with.
function [agents] = apply_comm_model_distance(SIM, agents, threshold)

    agents_arr = [];
    for id = 1:numel(agents)
        agent = agents{id};
        if SIM.OBJECTS(agent.objectID).type == 1
            agents_arr = [agents_arr, agent];
        end
    end
    for i = 1:numel(agents_arr)
        agents_arr(i).memory_id_comm = [];
        for j = 1:numel(agents_arr)
            if i ~= j
                detect = agents_arr(i).commRadius;
                Xi = SIM.OBJECTS(agents_arr(i).objectID).X;
                Xj = SIM.OBJECTS(agents_arr(j).objectID).X;
                if norm(Xj(1:3) - Xi(1:3)) <= detect
                    agents_arr(i).memory_id_comm = [agents_arr(i).memory_id_comm, j];
                end
            end
        end
    end
end

% Update the momory_id_obs list in each agent to include each agent they
% can observe.
function [agents] = apply_obs_model_distance(SIM, agents, threshold)

    agents_arr = [];
    for id = 1:numel(agents)
        agent = agents{id};
        if SIM.OBJECTS(agent.objectID).type == 1
            agents_arr = [agents_arr, agent];
        end
    end
    for i = 1:numel(agents_arr)
        agents_arr(i).memory_id_obs = [];
        for j = 1:numel(agents_arr)
            if i ~= j
                detect = agents_arr(i).obsRadius;
                Xi = SIM.OBJECTS(agents_arr(i).objectID).X;
                Xj = SIM.OBJECTS(agents_arr(j).objectID).X;
                if norm(Xj(1:3) - Xi(1:3)) <= detect
                    agents_arr(i).memory_id_obs = [agents_arr(i).memory_id_obs, j];
                end
            end
        end
    end
end

% Break agents up into groups based on communication graph
function [agent_groups] = break_agents_into_groups(SIM, agent_data)
    % Split into groups of agents that can communicate with eachother, 
    % and puts those agents in a group. Continues this along the chain
    % until there are no more agents in this group, then finds the other
    % isolated groups. 
    
    agent_groups = [];
    num_groups = 0;
    % While there are still agents to group up
    while numel(agent_data) > 0
        
        % Start with the first remaining agent
        group = agent_data{1};
        id_obs = group(1).memory_id_comm;
        new_group = 1;
        agent_data(1) = [];
        
        % while there are new agents in the group
        while new_group > 0
            
            % Get a list of the newly-observed IDs
            len = numel(group);
            tmp = group(1,len-new_group+1:len);
            id_obs_new = [];
            for m = 1:numel(tmp)
                id_obs_new = [id_obs_new, tmp(m).memory_id_comm];
            end
            id_obs_new = sort(unique(id_obs_new));
            id_obs_new = setdiff(id_obs_new, id_obs);
            id_obs = sort([id_obs, id_obs_new]);
            new_group = 0;
            indices = [];
            
            % Get the agents with ids matching the observed list
            for i = id_obs
                for j = 1:numel(agent_data)
                    if agent_data{j}.objectID == i
                        group = [group, agent_data{j}];
                        new_group = new_group + 1;
                        indices = [indices, j];
                    end
                end
            end
            
            % Remove grouped agents from the general pool
            for i = sort(indices, 'descend')
                agent_data(i) = [];
            end
        end
        agent_groups{num_groups+1} = group;
        num_groups = num_groups + 1;
    end
end

% Create a graph from a group of agents
function [graph, id_to_index] = create_graph(agents)
    adj = eye(numel(agents));
    
    id_to_index = [];
    for i = 1:numel(agents)
        id_to_index(i) = agents(i).objectID;
    end
    
    for i = 1:numel(agents)
        for j = agents(i).memory_id_comm
            adj(i,find(j == id_to_index)) = 1;
        end
    end
    
    if size(adj,1) == 1 & size(adj,2) == 1
        tmp = 0;
    end
    
    graph = generate_graph(adj);
    
end

% Perform one consensus step
function [consensus_data] = consensus_group(agents, step, num_steps)
    % Perform consensus computations
    
    %% Initialize Graph
    % Generate graph
        % Create adjacency matrix
        % use generate_graph function
    [graph, id_to_index] = create_graph(agents);
    size_comp = networkComponents(graph.p);
    
    %% Compute and store consensus variables
    for i = 1:numel(agents)
        
        % Grab variables from neighboring agents
        Y_local = [];
        y_local = [];
        idx_neighbors = agents(i).memory_id_comm;
        
        for j = 1:numel(idx_neighbors)
            agent = agent_with_id(agents, idx_neighbors(j));
            Y_local(:,:,j) = (agent.memory_Y);
            y_local(:,:,j) = (agent.memory_y);
        end
        
        % Compute and apply CI weights
        [weights_ci,Y_prior,y_prior] = calc_ci_weights_ver3(Y_local,y_local,'det');
        
        delta_I = zeros(size(agents(1).memory_I));
        delta_i = zeros(size(agents(1).memory_i));
        
        for j = 1:numel(agents)
            p_jk = graph.p(i,j);
            
            delta_I = delta_I + p_jk*agents(j).memory_I;
            delta_i = delta_i + p_jk*agents(j).memory_i;
        end
        
        ratio = step / num_steps;
        Y = Y_prior + ratio*size_comp(i)*delta_I;
        y = y_prior + ratio*size_comp(i)*delta_i;
        
        consensus_data{i}.Y = Y;
        consensus_data{i}.y = y;
        consensus_data{i}.Y_prior = Y_prior;
        consensus_data{i}.y_prior = y_prior;
        consensus_data{i}.delta_I = delta_I;
        consensus_data{i}.delta_i = delta_i;
    end
end

function [agent_groups] = consensus(agent_groups)
    num_steps = 20;
    for group_num = 1:numel(agent_groups)
        if numel(agent_groups{group_num}) == 1
            agent_groups{group_num}.memory_Y = agent_groups{group_num}.memory_Y + agent_groups{group_num}.memory_I;
            agent_groups{group_num}.memory_y = agent_groups{group_num}.memory_y + agent_groups{group_num}.memory_i;
            agent_groups{group_num}.memory_P = inv(agent_groups{group_num}.memory_Y);
            agent_groups{group_num}.memory_x = inv(agent_groups{group_num}.memory_Y) * agent_groups{group_num}.memory_y;
        else
            % Compute first consensus step
            step = 1;
            for i = 1:numel(agent_groups{group_num})
                consensus_data{step, group_num}{i}.Y_prior = agent_groups{group_num}(i).memory_Y;
                consensus_data{step, group_num}{i}.y_prior = agent_groups{group_num}(i).memory_y;
                consensus_data{step, group_num}{i}.delta_I = agent_groups{group_num}(i).memory_I;
                consensus_data{step, group_num}{i}.delta_i = agent_groups{group_num}(i).memory_i;
            end

            % Compute the remaining consensus steps
            for step = 2:num_steps
                consensus_data{step, group_num} = consensus_group(agent_groups{group_num}, step, num_steps);

                % After all agents' variables have been computed, store them
                for i = 1:numel(consensus_data{step, group_num})
                    agent_groups{group_num}(i).memory_Y = consensus_data{step, group_num}{i}.Y_prior;
                    agent_groups{group_num}(i).memory_y = consensus_data{step, group_num}{i}.y_prior;
                    agent_groups{group_num}(i).memory_I = consensus_data{step, group_num}{i}.delta_I;
                    agent_groups{group_num}(i).memory_i = consensus_data{step, group_num}{i}.delta_i;
                end
            end

            % Store final consensus in each agent
            for i = 1:numel(consensus_data{step, group_num})
                agent_groups{group_num}(i).memory_Y = consensus_data{step, group_num}{i}.Y;
                agent_groups{group_num}(i).memory_y = consensus_data{step, group_num}{i}.y;
                agent_groups{group_num}(i).memory_P = inv(consensus_data{step, group_num}{i}.Y);
                agent_groups{group_num}(i).memory_x = inv(consensus_data{step, group_num}{i}.Y) * consensus_data{step, group_num}{i}.y;
            end
        end
    end
end

% function [gt_estimation] = gt(agent_data)
%     y = y_1 + i_1 + i_2;
% end

function [rel_pos] = position_from_id(agent, id)

    x = agent.memory_x;
    index_id = find(id == agent.memory_id_list);
    index_agent = find(agent.objectID == agent.memory_id_list);
    dim_state = agent.dim_state;

    % Generate indices (to make the assignment setp shorter)
    meas_low = dim_state*(index_id - 1) + 1;
    meas_high = dim_state*index_id;
    agent_low = dim_state*(index_agent - 1) + 1;
    agent_high = dim_state*index_agent;

    rel_pos = x(meas_low:meas_high) - x(agent_low:agent_high);
    rel_pos = rel_pos(1:2);
end

% /////////////////// SIMULATION OUTPUT OPERATIONS ////////////////////////

% UPDATE THE OBJECT PROPERTIES
function [referenceObject,objectEVENTS] = UpdateObject(SIM,objectIndex,referenceObject)
% This function updates a referenceObject class against the rest of the
% objectIndex, independantly of the SIM.OBJECTS META data.
% INPUTS:
% SIM             - Current META structure
% objectIndex     - The current object list
% referenceObject - The agent object being updated

% OUTPUTS:
% referenceObject - The updated agent class

objectEVENTS = []; % Container for object based events

% GET THE OBJECTS EQUIVALENT SIM OBJECT
SIMfirstObject = SIM.OBJECTS(SIM.globalIDvector == referenceObject.objectID);
firstObject = objectIndex{SIM.globalIDvector == SIMfirstObject.objectID};

% If the observing object is "2D"
if referenceObject.Is3D()
    dimensionIndices = 1:3;
else
    dimensionIndices = 1:2;
end

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling ID:%d\t name: %s\n',SIMfirstObject.objectID,referenceObject.name);

end

% Default timing parameters
ENV = SIM.TIME;
ENV.outputPath = SIM.outputPath;

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling ID:%d\t name: %s\n',SIMfirstObject.objectID,referenceObject.name);
end

observationPacket = [];
% SWITCH BEHAVIOUR BASED ON OBJECT TYPE
switch SIMfirstObject.type
    case OMAS_objectType.agent
        % AGENT - SIMULATION/ENVIROMENTAL FEEDBACK REQUIRED %%%%%%%%%%%%%%
        % Container for agent detection packet
        
        % //// PARSE THE ENVIRONMENTAL OBJECTS OBSERVABLE TO THE AGENT ////
        for ID2 = firstObject.memory_id_obs
            % GET THE SIM OBJECT OF THE EVALUATION OBJECT
            SIMsecondObject = SIM.OBJECTS(ID2);
            
            % GET THE EQUIVALENT OBJECT TO EVALUATE AGAINST THE REFERENCE
            secondObject = objectIndex{SIM.globalIDvector == SIMsecondObject.objectID};
                            
            % ///////////////// ELSE AGENT IS DETECTED ////////////////////
            % ////// (GENERATE A LOCALISED PROXIMITY DESCIPTION) //////////
            % The second agent is within global seperation conditions. A
            % communication packet is therefore assembled containing the
            % data on the second object. rotated into the agent-local
            % coordinate frame.
            
            % GET THE RELATIVE POSITIONS OF THE AGENT AND OBJECT IN THE 
            % GLOBAL COORDINATE SYSTEM (x_rel = x_B - x_A)
            relativePosition = [position_from_id(firstObject, secondObject.objectID); 0];
            relativeVelocity = SIMsecondObject.globalState(4:6) - SIMfirstObject.globalState(4:6);
            % ROTATE THE GLOBAL STATE OF THE OBJECT INTO THE AGENT FRAME
            observedPosition = SIMfirstObject.R*relativePosition;          % Rotate the from the global into the body frame of the simReference
            observedVelocity = SIMfirstObject.R*relativeVelocity;
            % SPHERICAL REPRESENTATION
            observedRange     = norm(observedPosition);
            observedElevation = asin(observedPosition(3)/observedRange);
            observedHeading   = atan2(observedPosition(2),observedPosition(1));
            % PRESENT SIZE PROPERTIES
            observedRadius    = SIMfirstObject.radius;     
            observedAngularWidth = 2*asin(observedRadius/(observedRange + observedRadius));               
                        
            % OBJECT GEOMETRY PREPARATION
            % In order for the agent to observe the object correctly. The
            % geometry must be translated and rotated into the relative
            % frame of the detecting agent.
            % ASSUMPTION:
            % - The geometry is normalised to the body axes of the object.
            % TO DO:
            % - If the object is on the edge of vision.. only part of the
            %   object will be visable to the agent. a process must be in
            %   place to create a subset of the geometry.

            % THE GEOMETRIC PARAMETERS
            if size(secondObject.GEOMETRY.vertices,1) > 0
                % Rotate evaluation object geometry into the global space,
                % then rotate it back into the local space of the reference.                
                % THE RELATIVE ROTATIONS OF THE SECOND BODY
                relativeR = SIMfirstObject.R'*SIMsecondObject.R; 
                % THE CONTAINER FOR THE RELATIVE GEOMETRY
                observedGeometry = struct('vertices',secondObject.GEOMETRY.vertices*relativeR + observedPosition',...
                                          'normals', secondObject.GEOMETRY.normals*relativeR,...
                                          'faces',   secondObject.GEOMETRY.faces,...
                                          'centroid',secondObject.GEOMETRY.centroid + observedPosition');
                % [TO-DO] Reduce the structure sent to match the exposed geometry                      
                                      
%                 % PROCESS THE VISIBLE GEOMETRY (The sub-set of the geometry that is visible to the agent.)
%                 observedGeometry = OMAS_restrictedGeometry(zeros(3,1),SIMfirstObject.detectionRadius,relativeGeometry);                   
            else
                % THE OBSERVED GEOMETRY IS EMPTY
                observedGeometry = secondObject.GEOMETRY;                  % Pass the empty structure
            end
            
            % OBJECT PRIORITY (IF WAYPOINT)
            observedPriority = NaN;
            if SIMsecondObject.type == OMAS_objectType.waypoint
                association = objectIndex{ID2}.GetAgentAssociation(referenceObject);
                observedPriority = association.priority;
            end
            
            % A DEBUG STUCTURE 
            DEBUG = struct('globalPosition',SIMsecondObject.globalState(dimensionIndices),...     % Added for simplicity
                           'globalVelocity',SIMsecondObject.globalState(dimensionIndices + 3),... % Added for simplicity
                           'priority',observedPriority);                                          % The objects global priority for that agent
            
            % ///////////////// PREPARE DETECTION PACKET //////////////////
            detectionObject = struct('objectID',SIMsecondObject.objectID,...            % The object ID (observed)
                                     'name',SIMsecondObject.name,...                    % The object name tag (observed)
                                     'type',SIMsecondObject.type,...                    % The objects sim-type enum
                                     'radius',observedRadius,...                        % The objects true size
                                     'position',observedPosition(dimensionIndices,1),...% The apparent position in the relative frame
                                     'velocity',observedVelocity(dimensionIndices,1),...% The apparent velocity in the relative frame
                                     'range',observedRange,...                          % The apparent range
                                     'elevation',observedElevation,...                  % The apparent inclination angle
                                     'heading',observedHeading,...                      % The apparent Azimuth angle
                                     'width',observedAngularWidth,...                   % The apparent angular width at that range    
                                     'geometry',observedGeometry,...                    % The observable geometrical components
                                     'colour',SIMsecondObject.colour,...                % Finally the simulation's colourID 
                                     'DEBUG',DEBUG);                                    % Pass additional parameters (not otherwise known)
            observationPacket = vertcat(observationPacket,detectionObject);             % Append object to packet to agent
        end

        % /////////////// COMPUTE THE LOCAL AGENT CYCLE ///////////////////
        % Given the detection object defined for this agent, compute the
        % agents cycle with this information.
        % SEND OBJECT OBSERVERATION PACKET TO AGENT
        referenceObject.main(ENV,observationPacket);
    case OMAS_objectType.obstacle
        % PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        referenceObject = referenceObject.main(ENV);
    case OMAS_objectType.waypoint
        % WAYPOINTS ARE CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%
        referenceObject = referenceObject.main(ENV);
    otherwise
        % CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%
        % DO NOT UPDATE
end
end

% /////////////////// SIMULATION OUTPUT OPERATIONS ////////////////////////
% GET INITIAL OUTPUT DATA STRUCTURE
function [DATA] = GetOutputStructure(SIM)
% This function assembles the initial output DATA structure from the
% simulation input variables and it back as a simulation output container.
% INPUT:
% SIM  - A local reference to the META structure
% OUPUT:
% DATA - The output data structure

% CALCULATE THE INITIAL DATA PARAMETERS
globalStateNum = size(SIM.OBJECTS(1).globalState,1);
systemStates = globalStateNum*SIM.totalObjects;                            % The total number of system states
% DEFINE THE GLOBAL STATE INDICES
id_initial = zeros(SIM.totalAgents, SIM.totalObjects, SIM.TIME.numSteps);
comm_initial = zeros(SIM.totalAgents, SIM.totalAgents, SIM.TIME.numSteps);
estimate_initial = zeros(SIM.totalAgents, SIM.totalObjects, 12, SIM.TIME.numSteps);
covariance_initial = zeros(SIM.totalAgents, SIM.totalObjects, 12, 12, SIM.TIME.numSteps);
indexSet = zeros(2,length(SIM.OBJECTS));                                   % Indices are defined [start;end]*n
indexSet(1,:) = (0:globalStateNum:systemStates-globalStateNum) + 1;        % Declare the state indices
indexSet(2,:) = globalStateNum:globalStateNum:systemStates;

% GENERATE THE OUTPUT STRUCTURE FROM THE SIMULATION INPUT PARAMETERS
DATA = struct('outputPath',[SIM.outputPath,'DATA.mat'],...
               'outputDir',SIM.outputPath,...
            'totalObjects',SIM.totalObjects,...
             'totalAgents',SIM.totalAgents,...
          'totalObstacles',SIM.totalObstacles,...
          'totalWaypoints',SIM.totalWaypoints,...
              'timeVector',SIM.TIME.timeVector,...
                      'dt',SIM.TIME.dt,...
              'stateIndex',indexSet,...
                     'ids',id_initial,...
            'Observations',id_initial,...
             'Connections',comm_initial,...
           'estimates_rel',estimate_initial,...
     'estimate_errors_rel',estimate_initial,...
               'estimates',estimate_initial,...
         'estimate_errors',covariance_initial,...
     'estimate_covariance',[],...
      'globalTrajectories',NaN(systemStates,SIM.TIME.numSteps));         % Prepare the output container
%       'globalTrajectories',zeros(systemStates,SIM.TIME.numSteps));         % Prepare the output container
end