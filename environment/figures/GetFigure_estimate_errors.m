function [currentFigure,figureHandle] = GetFigure_estimate_errors(SIM,objectIndex,DATA,currentFigure)

% CONFIGURE THE PLOT ATTRIBUTES
figurePath = strcat(SIM.outputPath,'isometric');
figureHandle = figure('Name','OpenMAS estimate errors view');
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW
ax = axes(figureHandle);
legendCounter = 1; legendEntries = cell(DATA.totalAgents^2 - DATA.totalAgents,1);

hold on; grid on; box on; grid minor;

for ID1 = 1:DATA.totalAgents 
    idleFlag = NaN('double');
    [objectStates] = OMAS_getTrajectoryData(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,idleFlag);
    agent_pos = objectStates(1:3,:);
    
    % Poorly written code to dynamically find the start and end of the
    % disconnection of a 3-agent system. 
    connected_time = 1 - squeeze(DATA.Connections(2, 2, :) / 3);
    connected_time = connected_time'.*(1:numel(connected_time));
    
    if sum(connected_time) > 0
        start = min(connected_time(connected_time > 0)) * DATA.dt;
        finish = max(connected_time) * DATA.dt;
        
        y_max = 1.2;
        x = [start, start, finish, finish];
        y = [0, y_max, y_max, 0];
        fill(x, y, 'k', 'facealpha', 0.2,'edgecolor','none','HandleVisibility','off');
    end
    
    for ID2 = 1:DATA.totalObjects
        if SIM.OBJECTS(ID2).type == 1 && ID1 ~= ID2
            % GET OBJECT OVERVIEW DATA
            legendString = sprintf('[Agent-%s] estimating [Agent-%s]',num2str(SIM.OBJECTS(ID1).objectID),num2str(SIM.OBJECTS(ID2).objectID));
            legendEntries(legendCounter) = {legendString};

            positions = squeeze(DATA.estimate_errors(ID1, ID2, 1:3, :));
            time = (1:size(positions,2)) * DATA.dt;
            sim_time = 1:(SIM.TIME.endTime / DATA.dt);
            % EXTRACT STATE TIME-SERIES DATA UPTO THE IDLE POINT
            plot(time(sim_time), (positions(1,sim_time).^2 + positions(2,sim_time).^2).^0.5,...
                  'LineWidth',DATA.figureProperties.lineWidth,...
                  'Color',SIM.OBJECTS(ID1).colour);

            legendCounter = legendCounter + 1;
        end
    end
end

legend('off')
hold on;

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

% xlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
% ylim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
%zlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
%xlim(ax,[DATA.figureProperties.axisMinimums(1)-0.1,DATA.figureProperties.axisMaximums(1)+0.1]);
%ylim(ax,[DATA.figureProperties.axisMinimums(2)-0.1,DATA.figureProperties.axisMaximums(2)+0.1]);
% zlim(ax,[DATA.figureProperties.axisMinimums(3),DATA.figureProperties.axisMaximums(3)]);

%axis vis3d equal;     
%view([-24 36]);
%set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
grid on; grid minor;
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath); 

% Publish as .pdf if requested
if DATA.figureProperties.publish
	GetFigurePDF(figureHandle,figurePath);   
end

% FIGURE COMPLETE
currentFigure = currentFigure + 1;
% CLEAN UP
clearvars -except currentFigure figureHandle
end
