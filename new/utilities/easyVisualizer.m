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
        function fig = plotFillElementData(obj, elementcoordinates, elementfillstatus, results, timestep, colorlimit, filled, crossSection, displayfig, crossSectionValue)
            [filledelements, emptyelements] = obj.getFilledElements(elementfillstatus);
        
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);
        
            filledIdx = filledelements{timestep};
            emptyIdx = emptyelements{timestep};
            if nargin < 7
                filled = true;
            end
            if nargin < 8
                crossSection = '';
            end
            
            if nargin < 9
                displayfig = false;
            end

            if nargin < 10
                switch crossSection
                    case 'XY'
                        crossSectionValue = (max(z) + min(z)) / 2;
                    case 'XZ'
                        crossSectionValue = (max(y) + min(y)) / 2;
                    case 'YZ'
                        crossSectionValue = (max(x) + min(x)) / 2;
                    otherwise 
                        crossSectionValue = 0;
                end
        
            end
        
            % Apply cross-section filter if specified
            if ~isempty(crossSection)
                if displayfig == true
                    fig = figure('Units', 'normalized', 'OuterPosition', [0.1 0.1 0.3 0.8]);
                else
                    fig = figure('Units', 'normalized', 'OuterPosition', [0.1 0.1 0.3 0.8], 'Visible', 'off');
                end
                tiledlayout(2, 1);
                switch crossSection
                    case 'XY'
                        filterIdx = z <= crossSectionValue;
                    case 'XZ'
                        filterIdx = y >= crossSectionValue;
                    case 'YZ'
                        filterIdx = x >= crossSectionValue;
                    otherwise
                        error('Invalid cross-section specified.');
                end
                filledIdx = filledIdx(filterIdx(filledIdx));
                emptyIdx = emptyIdx(filterIdx(emptyIdx));
            else
                if displayfig == true
                    fig = figure();
                else
                    fig = figure('Visible', 'off')
                end
            end
            
            if ~isempty(crossSection)
                nexttile;
            end
        
            if filled
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
            
            % Plot 3D cross-sectional view if specified
            if ~isempty(crossSection)
                nexttile;
                scatter3(x(filledIdx), y(filledIdx), z(filledIdx), 20, results(filledIdx, timestep), 'filled', 'd', 'MarkerEdgeColor', 'none');
                hold on
                scatter3(x(emptyIdx), y(emptyIdx), z(emptyIdx), 20, 'd', 'w', 'MarkerFaceAlpha', 0.01,'MarkerEdgeAlpha', 0.01);
        
                colormap('jet')
                if strcmp(colorlimit, 'yes')
                    clim([min(results(:)), max(results(:))])
                end
                colorbar
                axis equal
                hold off
            end
        end

        % Function for plotting elements to refine and elements that remain
        % unchanged side by side
        function fig = plotElementsToRefine(obj, elementcoordinates, highdiffelements, normalelements, affectedelements, crossSection, displayfig, crossSectionValue)
            x = elementcoordinates(:, 1);
            y = elementcoordinates(:, 2);
            z = elementcoordinates(:, 3);
            
            if nargin < 6
                crossSection = '';
            end

            if nargin < 7
                displayfig = false;
            end
        
            if nargin < 8
                switch crossSection
                    case 'XY'
                        crossSectionValue = (max(z) + min(z)) / 2;
                    case 'XZ'
                        crossSectionValue = (max(y) + min(y)) / 2;
                    case 'YZ'
                        crossSectionValue = (max(x) + min(x)) / 2;
                    otherwise 
                        crossSectionValue = 0;
                end
            end
            
            % Apply cross-section filter if specified
            if ~isempty(crossSection)
                if displayfig == true
                    fig = figure('Units', 'normalized', 'OuterPosition', [0.1 0.1 0.3 0.8]);
                else
                    fig = figure('Units', 'normalized', 'OuterPosition', [0.1 0.1 0.3 0.8], 'Visible','off');
                end
                tiledlayout(2, 1);
                switch crossSection
                    case 'XY'
                        filterIdx = z <= crossSectionValue;
                    case 'XZ'
                        filterIdx = y >= crossSectionValue;
                    case 'YZ'
                        filterIdx = x >= crossSectionValue;
                    otherwise
                        error('Invalid cross-section specified.');
                end
                highdiffelements_cs = highdiffelements(filterIdx(highdiffelements(:, 1)), :);
            else
                if displayfig == true
                    fig = figure();
                else 
                    fig = figure('Visible','off');
                end
            end
        
        
            if ~isempty(crossSection)
                nexttile;
            end
            
            scatter3(x(highdiffelements(:, 1)), y(highdiffelements(:, 1)), z(highdiffelements(:, 1)), 100, 'filled', "d", "red");
            hold on
        
            if affectedelements ~= 0
                affectedelements_disp = affectedelements(find(~ismember(affectedelements, highdiffelements(:, 1))));
                unchangedelements_disp = normalelements(find(~ismember(normalelements, affectedelements)));
                scatter3(x(affectedelements_disp), y(affectedelements_disp), z(affectedelements_disp), 100, 'filled', 'd', 'green', 'MarkerFaceAlpha', 0.5,'MarkerEdgeAlpha', 0.8);
                if ~isempty(crossSection)
                    affectedelements_disp_cs = affectedelements_disp(filterIdx(affectedelements_disp));
                end
            else
                unchangedelements_disp = normalelements;
            end
        
            if ~isempty(crossSection)
                unchangedelements_disp_cs = unchangedelements_disp(filterIdx(unchangedelements_disp));
            end
        
            scatter3(x(unchangedelements_disp), y(unchangedelements_disp), z(unchangedelements_disp), 100,'filled', "d","y", 'MarkerFaceAlpha', 0.15,'MarkerEdgeAlpha', 0.2);
            axis equal
            hold off
        
            if ~isempty(crossSection)
                nexttile;
                scatter3(x(highdiffelements_cs(:, 1)), y(highdiffelements_cs(:, 1)), z(highdiffelements_cs(:, 1)), 100, 'filled', "d", "red");
                hold on
                scatter3(x(unchangedelements_disp_cs), y(unchangedelements_disp_cs), z(unchangedelements_disp_cs), 100,'filled', "d","y", 'MarkerFaceAlpha', 0.15,'MarkerEdgeAlpha', 0.2);
                axis equal
                hold off
            end
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