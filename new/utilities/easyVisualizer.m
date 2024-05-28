classdef easyVisualizer
    % class containing functions for easy visualization of results and other
    % data
    methods
    % Below here are functions for handy data visualization
    % Function for plotting results

            % Function for separating the indices of filled and empty elements for
        % every timestep
        function [filledelements, emptyelements] = getFilledElements(obj, elementfillstatus)
            filledelements = cell(1, size(elementfillstatus, 2));
            emptyelements = cell(1, size(elementfillstatus, 2));
            
            for i = 1:size(elementfillstatus, 2)
                filledindices = find(elementfillstatus(:, i) >= 1);
                emptyindices = find(elementfillstatus(:, i) < 1);
                filledelements{i} = filledindices;
                emptyelements{i} = emptyindices;
            end
        end

        % Function for calculating the postions of the elements (their
        % centroids)
        function [elementcoordinates] = getElementCoordinates(obj, elementnodeids, nodepos)
            elementcoordinates = zeros(size(elementnodeids, 1), 3);
            for i = 1:size(elementnodeids, 1)
                columnids = elementnodeids(i, :);
                x = mean(nodepos(columnids, 1));
                y = mean(nodepos(columnids, 2));
                z = mean(nodepos(columnids, 3));
                elementcoordinates(i, :) = [x, y, z];
            end
        end

        function plotElementData(obj, elementcoordinates, results, timestep, colorlimit)
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);
        
            scatter3(x, y, z, 20, results(:, timestep), "filled", 'd');
            colormap('jet')
            axis equal
            if strcmp(colorlimit, 'yes')
                clim([min(results(:)), max(results(:))])
            end
            colorbar
        end


        % Function for plotting the filling method while also displaying the 
        % element values for a feature in the filled region. The flowfront is
        % considered as empty so element values in the flowfront are ignored
        function plotFillElementData(obj, elementcoordinates, elementfillstatus, results, timestep, colorlimit, filled)
            if nargin < 7
                filled = true;
            end
            [filledelements, emptyelements] = obj.getFilledElements(elementfillstatus);
        
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);
            
            if filled == true
                scatter3(x(filledelements{timestep}), y(filledelements{timestep}), z(filledelements{timestep}), 20, results(filledelements{timestep}, timestep), 'filled', 'd', 'MarkerEdgeColor', 'none');
            else
                scatter3(x(filledelements{timestep}), y(filledelements{timestep}), z(filledelements{timestep}), 20, results(filledelements{timestep}, timestep),'d');
            end
            colormap('jet')
            if strcmp(colorlimit, 'yes')
                clim([min(results(:)), max(results(:))])
            end
            colorbar
            hold on
            scatter3(x(emptyelements{timestep}), y(emptyelements{timestep}), z(emptyelements{timestep}), 20, 'd', 'w', 'MarkerFaceAlpha', 0.01,'MarkerEdgeAlpha', 0.01);
            axis equal
            hold off
        end

        % Function for plotting elements to refine and elements that remain
        % unchanged side by side
        function fig = plotElementsToRefine(obj, elementcoordinates, highdiffelements, normalelements, affectedelements)
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);

            % if nargout == 0
            %     fig = figure;
            % else
            %     fig = figure;
            % end

            fig = figure;

            scatter3(x(highdiffelements(:, 1)), y(highdiffelements(:, 1)), z(highdiffelements(:, 1)), 100, 'filled', "d", "red", 'MarkerFaceAlpha', 0.8,'MarkerEdgeAlpha', 0.8);
            hold on

            if affectedelements ~= 0
                affectedelements_disp = affectedelements(find(~ismember(affectedelements, highdiffelements(:, 1))));
                unchangedelements_disp = normalelements(find(~ismember(normalelements, affectedelements)));
                scatter3(x(affectedelements_disp), y(affectedelements_disp), z(affectedelements_disp), 100, 'filled', 'd', 'green', 'MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha', 0.8);
            else
                unchangedelements_disp = normalelements;
            end
            scatter3(x(unchangedelements_disp), y(unchangedelements_disp), z(unchangedelements_disp), 100,'filled', "d","y", 'MarkerFaceAlpha', 0.1,'MarkerEdgeAlpha', 0.2);
            axis equal
            hold off
        end

        % Function for plotting new elements and old elements (that
        % remained unchanged) side by side
        function plotNewElements(obj, elementnodeids, newelementnodeids, affectedelements, newelementcoordinates)
            unchangedelements = elementnodeids;
            unchangedelements(affectedelements, :) = [];
            oldelements = find(ismember(newelementnodeids, unchangedelements, 'rows'));
            newelements = find(~ismember(newelementnodeids, unchangedelements, 'rows'));
    
            x = newelementcoordinates(:, 1);
            y = newelementcoordinates(:, 2);
            z = newelementcoordinates(:, 3);

            scatter3(x(oldelements), y(oldelements), z(oldelements), 100,  "filled", "d","y", 'MarkerFaceAlpha', 0.1,'MarkerEdgeAlpha', 0.2);
            hold on
            scatter3(x(newelements), y(newelements), z(newelements), 100, 'filled', 'd', 'green');
            axis equal
            hold off
        end
    end
end