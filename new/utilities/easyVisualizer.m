classdef easyVisualizer
    % class containing functions for easy visualization of results and other
    % data
    methods
    % Below here are functions for handy data visualization
    % Function for plotting results
        function plotElementData(elementcoordinates, results, timestep, colorlimit)
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);
        
            scatter3(x, y, z, 20, results(:, timestep), "filled")
            axis equal
            if strcmp(colorlimit, 'yes')
                clim([min(results(:)), max(results(:))])
            end
            colorbar
        end


        % Function for plotting the filling method while also displaying the 
        % element values for a feature in the filled region. The flowfront is
        % considered as empty so element values in the flowfront are ignored
        function plotFillElementData(elementcoordinates, elementfillstatus, results, timestep, colorlimit)
            [filledelements, emptyelements] = getFilledElements(elementfillstatus);
        
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);
        
            scatter3(x(filledelements{timestep}), y(filledelements{timestep}), z(filledelements{timestep}), 20, results(filledelements{timestep}, timestep), "filled", 'MarkerEdgeColor', 'none')
            if strcmp(colorlimit, 'yes')
                clim([min(results(:)), max(results(:))])
            end
            colorbar
            hold on
            scatter3(x(emptyelements{timestep}), y(emptyelements{timestep}), z(emptyelements{timestep}), 20, 'filled', 'w', 'MarkerFaceAlpha', 0.2,'MarkerEdgeAlpha', 0.2)
            axis equal
            hold off
        end
    end
end