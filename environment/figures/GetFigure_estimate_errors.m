function [currentFigure,figureHandle] = GetFigure_estimate_errors(SIM,objectIndex,DATA,currentFigure)

% % CONFIGURE THE PLOT ATTRIBUTES
% figurePath = strcat(SIM.outputPath,'isometric');
% figureHandle = figure('Name','OpenMAS estimate errors view');
% set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
% set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% % setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW
% ax = axes(figureHandle);
% legendCounter = 1; legendEntries = cell(DATA.totalAgents^2 - DATA.totalAgents,1);
% 
% hold on; grid on; box on; grid minor;
% 
% for ID1 = 1:DATA.totalAgents 
%     
% %     % Poorly written code to dynamically find the start and end of the
% %     % disconnection of a 3-agent system. 
% %     connected_time = 1 - squeeze(DATA.Connections(2, 2, :) / 3);
% %     connected_time = connected_time'.*(1:numel(connected_time));
% %     
% %     if sum(connected_time) > 0
% %         start = min(connected_time(connected_time > 0)) * DATA.dt;
% %         finish = max(connected_time) * DATA.dt;
% %         
% %         y_max = 1.2;
% %         x = [start, start, finish, finish];
% %         y = [0, y_max, y_max, 0];
% %         fill(x, y, 'k', 'facealpha', 0.2,'edgecolor','none','HandleVisibility','off');
% %     end
%     
%     for ID2 = 1:DATA.totalObjects
%         if SIM.OBJECTS(ID2).type == 1 && ID1 ~= ID2
%             % GET OBJECT OVERVIEW DATA
%             legendString = sprintf('[Agent-%s] estimating [Agent-%s]',num2str(SIM.OBJECTS(ID1).objectID),num2str(SIM.OBJECTS(ID2).objectID));
%             legendEntries(legendCounter) = {legendString};
% 
%             positions = squeeze(DATA.estimate_errors(ID1, ID2, 1:3, :));
%             time = (1:size(positions,2)) * DATA.dt;
%             sim_time = 1:(SIM.TIME.endTime / DATA.dt);
%             errors = (positions(1,sim_time).^2 + positions(2,sim_time).^2).^0.5;
%             % EXTRACT STATE TIME-SERIES DATA UPTO THE IDLE POINT
%             plot(time(sim_time), errors,...
%                   'LineWidth',DATA.figureProperties.lineWidth,...
%                   'Color',SIM.OBJECTS(ID1).colour);
% 
%             legendCounter = legendCounter + 1;
%         end
%     end
% end
% 
% legend('off')
% hold on;
% 
% % Title
% title(ax,...
%     sprintf('Estimate errors over a period of %ss',num2str(SIM.TIME.endTime)),...
%     'interpreter',DATA.figureProperties.interpreter,...
%     'fontname',DATA.figureProperties.fontName,...
%     'fontweight',DATA.figureProperties.fontWeight,...
%     'fontsize',DATA.figureProperties.titleFontSize,...
%     'fontsmoothing','on');
% % X-Label
% xlabel(ax,'time (seconds)',...
%     'Interpreter',DATA.figureProperties.interpreter,...
%     'fontname',DATA.figureProperties.fontName,...
%     'Fontweight',DATA.figureProperties.fontWeight,...
%     'FontSize',DATA.figureProperties.axisFontSize,...
%     'FontSmoothing','on');
% % Y-Label
% ylabel(ax,'position estimation error magnitude ||E|| (m)',...
%     'Interpreter',DATA.figureProperties.interpreter,...
%     'fontname',DATA.figureProperties.fontName,...
%     'Fontweight',DATA.figureProperties.fontWeight,...
%     'FontSize',DATA.figureProperties.axisFontSize,...
%     'FontSmoothing','on');
% % Axes
% set(ax,...
%     'TickLabelInterpreter',DATA.figureProperties.interpreter,...
%     'fontname',DATA.figureProperties.fontName,...
%     'Fontweight',DATA.figureProperties.fontWeight,...
%     'FontSize',DATA.figureProperties.axisFontSize,...
%     'FontSmoothing','on',...
%     'Color',DATA.figureProperties.axesColor,...
%     'GridLineStyle','--',...
%     'GridAlpha',0.25,...
%     'GridColor','k');
% % Legend
% legend(legendEntries,...
%     'location','northeastoutside',...
%     'fontname',DATA.figureProperties.fontName,...
%     'interpreter',DATA.figureProperties.interpreter);

% xlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
% ylim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
%zlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
%xlim(ax,[DATA.figureProperties.axisMinimums(1)-0.1,DATA.figureProperties.axisMaximums(1)+0.1]);
%ylim(ax,[DATA.figureProperties.axisMinimums(2)-0.1,DATA.figureProperties.axisMaximums(2)+0.1]);
% zlim(ax,[DATA.figureProperties.axisMinimums(3),DATA.figureProperties.axisMaximums(3)]);

%axis vis3d equal;     
%view([-24 36]);
%set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
% grid on; grid minor;
% hold off;
% 
% % SAVE THE OUTPUT FIGURE
% savefig(figureHandle,figurePath); 
% 
% % Publish as .pdf if requested
% if DATA.figureProperties.publish
% 	GetFigurePDF(figureHandle,figurePath);   
% end
% 
% % FIGURE COMPLETE
% currentFigure = currentFigure + 1;

set_maxmin = 0;
for ID1 = 1:DATA.totalAgents 
    
    for ID2 = 1:DATA.totalAgents
        if SIM.OBJECTS(ID2).type == 1 && ID1 ~= ID2
            positions = squeeze(DATA.estimate_errors(ID1, ID2, 1:3, :));
            sim_time = 1:(SIM.TIME.endTime / DATA.dt);
            errors = (positions(1,sim_time).^2 + positions(2,sim_time).^2).^0.5;
            
            if ~set_maxmin
                err_max = max(errors,[],2)';
                set_maxmin = 1;
            else
                err_max = max([err_max; max(errors,[],2)']);
            end
        end
    end
end

err_max = 5*ceil(err_max/5);

for ID1 = 1:DATA.totalAgents 
    
    figureHandle = figure('Name',sprintf('OpenMAS agent %d estimated trajectories view', ID1));
    set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
    set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
    % setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW
    ax = axes(figureHandle);
    legendCounter = 1; legendEntries = cell(DATA.totalAgents - 1,1);

    hold on; grid on; box on; grid minor;
    
%     % THE ASSOCIATED LOGIC
%     objectID1 = objectIndex{SIM.globalIDvector == SIM.OBJECTS(ID1).objectID};
%     % EXTRACT FINAL POSITION DATA FROM THE TRAJECTORY MATRIX
%     [finalStates] = OMAS_getTrajectoryData(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,SIM.TIME.endStep);
%     finalPosition = finalStates(1:3,:);
%     
%     % LOCAL FIXED TO GLOBAL ROTATED
%     R_final = OMAS_geometry.quaternionToRotationMatrix(finalStates(7:10));
%     
%     % DISPLAY THE OBJECT
%     if numel(objectID1.GEOMETRY.vertices) > 0
%         patch(ax,'Vertices',objectID1.GEOMETRY.vertices*R_final + finalPosition',...
%             'Faces',objectID1.GEOMETRY.faces,...
%             'FaceColor',SIM.OBJECTS(ID1).colour,...
%             'EdgeColor',DATA.figureProperties.edgeColor,...
%             'EdgeAlpha',DATA.figureProperties.edgeAlpha,...  
%             'FaceLighting',DATA.figureProperties.faceLighting,...
%             'FaceAlpha',DATA.figureProperties.faceAlpha,...
%             'LineWidth',DATA.figureProperties.patchLineWidth);             % Patch properties            % Patch properties
%     else
%         % PLOT THE TERMINAL POSITIONS
%         plot3(ax,finalPosition(1),finalPosition(2),finalPosition(3),...
%               'Marker',SIM.OBJECTS(ID1).symbol,...
%               'MarkerSize',DATA.figureProperties.markerSize,...
%               'MarkerFaceColor',SIM.OBJECTS(ID1).colour,...
%               'MarkerEdgeColor',DATA.figureProperties.markerEdgeColor,...
%               'LineWidth',DATA.figureProperties.lineWidth,...
%               'LineStyle',DATA.figureProperties.lineStyle,...
%               'Color',SIM.OBJECTS(ID1).colour); 
%     end  
    
%     legendString = sprintf('[Agent-%s] trajectory',num2str(SIM.OBJECTS(ID1).objectID));
%     legendEntries(legendCounter) = {legendString};
%     
%     idleFlag = NaN('double');
%     [objectStates] = OMAS_getTrajectoryData(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,idleFlag);
%     agent_pos = objectStates(1:3,:);
% 
%     positions = agent_pos + squeeze(DATA.estimates_rel(ID1, ID1, 1:3, :));
%     time = (1:size(positions,2)) * DATA.dt;
%     sim_time = 1:(SIM.TIME.endTime / DATA.dt);
%     % EXTRACT STATE TIME-SERIES DATA UPTO THE IDLE POINT
%     plot(positions(1,:),positions(2,:),...
%           'LineStyle',DATA.figureProperties.lineStyle,...
%           'LineWidth',DATA.figureProperties.lineWidth,...
%           'Color',SIM.OBJECTS(ID1).colour);
% 
%     legendCounter = legendCounter + 1;
    
    for ID2 = 1:DATA.totalObjects
        if SIM.OBJECTS(ID2).type == 1 && ID1 ~= ID2
            % GET OBJECT OVERVIEW DATA
            legendString = sprintf('[Agent-%s] estimating [Agent-%s]',num2str(SIM.OBJECTS(ID1).objectID),num2str(SIM.OBJECTS(ID2).objectID));
            legendEntries(legendCounter) = {legendString};

            positions = squeeze(DATA.estimate_errors(ID1, ID2, 1:3, :));
            time = (1:size(positions,2)) * DATA.dt;
            sim_time = 1:(SIM.TIME.endTime / DATA.dt);
            errors = (positions(1,sim_time).^2 + positions(2,sim_time).^2).^0.5;
            cov = squeeze(DATA.estimate_covariance(ID2, ID1, 1:2, 1:2, sim_time));
            cov_mag = squeeze((cov(1,1,:).^2 + cov(1,2,:).*cov(2,1,:) + cov(2,2,:).^2).^0.5);
            % EXTRACT STATE TIME-SERIES DATA UPTO THE IDLE POINT
            d = errorbar(time(sim_time), errors, cov_mag,...
                  'LineStyle',DATA.figureProperties.lineStyle,...
                  'LineWidth',DATA.figureProperties.lineWidth,...
                  'Color',SIM.OBJECTS(ID2).colour);
            d.Bar.LineStyle = 'dotted';

            legendCounter = legendCounter + 1;
            
%             for t = sim_time
%                 if mod(t, round(5/DATA.dt)) == 0
%                     mu = positions(1:2,t);
%                     cov = squeeze(DATA.estimate_covariance(ID2, ID1, 1:2, 1:2, t));
%                     ellipse = error_ellipse(cov, mu);
%                 end
%             end
        end
    end
    
    legend('off')
%     hold on;

    % Title
    title(ax,...
        sprintf('Estimate errors over a period of %ss',num2str(SIM.TIME.endTime)),...
        'interpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
        'fontweight',DATA.figureProperties.fontWeight,...
        'fontsize',DATA.figureProperties.titleFontSize,...
        'fontsmoothing','on');
    % X-Label
    xlabel(ax,'time (seconds)',...
        'Interpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on');
    % Y-Label
    ylabel(ax,'position estimation error magnitude ||E|| (m)',...
        'Interpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on');
    % Axes
    set(ax,...
        'TickLabelInterpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on',...
        'Color',DATA.figureProperties.axesColor,...
        'GridLineStyle','--',...
        'GridAlpha',0.25,...
        'GridColor','k');
    % Legend
    legend(legendEntries,...
        'location','northeastoutside',...
        'fontname',DATA.figureProperties.fontName,...
        'interpreter',DATA.figureProperties.interpreter);

    ylim(ax,[0,err_max]);
    %zlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
    %xlim(ax,[DATA.figureProperties.axisMinimums(1)-0.1,DATA.figureProperties.axisMaximums(1)+0.1]);
    %ylim(ax,[DATA.figureProperties.axisMinimums(2)-0.1,DATA.figureProperties.axisMaximums(2)+0.1]);
    % zlim(ax,[DATA.figureProperties.axisMinimums(3),DATA.figureProperties.axisMaximums(3)]);

    % axis vis3d equal;     
    % view([-24 36]);
    % set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
    grid on; grid minor;
    hold off;
end

% CLEAN UP
clearvars -except currentFigure figureHandle
end
