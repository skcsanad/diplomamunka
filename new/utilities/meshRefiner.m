classdef meshRefiner
    % class containing functions for doing the actual mesh refinement

    methods
        function elementvalues = getElementValues(obj, elementnodeids, values) % Function for mapping node values to elements
            elementvalues = zeros(size(elementnodeids, 1), size(elementnodeids, 2), size(values, 1));
            for i = 1:size(elementnodeids, 1)
                columnids = elementnodeids(i, :);
                currentvalues = values(:, columnids);
                elementvalues(i, :, :) = currentvalues';
            end
        end


        % Function for calculating de diffs in the elements
        function elementdiffs = calcElementDiffs(obj, elementvalues, mode)
            % Default mode: range
            if nargin < 3
                mode = 'range';
            end
            % Different diffmetrics
            if strcmp(mode, 'range')
                elementdiffs = (max(elementvalues, [], 2) - min(elementvalues, [], 2)); %table2array(rowfun(@(x) max(x) - min(x), table(elementvalues)));
            elseif strcmp(mode, 'normrange')
                elementdiffs = (max(elementvalues, [], 2) - min(elementvalues, [], 2)) ./ abs(mean(elementvalues, 2));
            elseif strcmp(mode, 'relrange')
                elementdiffs = (max(elementvalues, [], 2) - min(elementvalues, [], 2)) ./ max(elementvalues, [], 2);
            end
            elementdiffs = reshape(elementdiffs, size(elementdiffs, 1), size(elementdiffs, 3));
            elementdiffs(isnan(elementdiffs)) = 0;
        end
        
        
        % Function for calculating the fillstatus in the elements possible node
        % values of 3 and 4 are ignored because they are not present in the data
        % and will likely not be used anyways
        function elementfillstatus = calcElementFillstatus(obj, elementvalues)
            elementfillstatus = zeros(size(elementvalues, 1), size(elementvalues, 3));
            for i = 1:size(elementvalues, 1)
                for j = 1:size(elementvalues, 3)
                    if all(elementvalues(i, :, j) == 0)
                        elementfillstatus(i, j) = 1;
                    else
                        elementfillstatus(i, j) = 0;
                    end
                end
            end
        end
        
        
        % Function for separating the indices of filled and empty elements for
        % every timestep
        function [filledelements, newfilledelements, emptyelements] = getFilledElements(obj, elementfillstatus)
            filledelements = cell(1, size(elementfillstatus, 2));
            emptyelements = cell(1, size(elementfillstatus, 2));
            allfilledelements = [];
            newfilledelements = cell(1, size(elementfillstatus, 2));
            
            for i = 1:size(elementfillstatus, 2)
                filledindices = find(elementfillstatus(:, i) >= 1);
                newfilledindices = filledindices(~ismember(filledindices, allfilledelements));
                allfilledelements = [allfilledelements; filledindices];
                emptyindices = find(elementfillstatus(:, i) < 1);
                filledelements{i} = filledindices;
                newfilledelements{i} = newfilledindices;
                emptyelements{i} = emptyindices;
            end
        end
        
        
        % Function for retrieving elements with high diff values, the normal
        % elements and a cell containing the highdiff elements by timestep - this
        % might be removed if there is no need for use
        function [highdiffelements, normalelements, highdiffsbytimestep] = calcHighDiffElements(obj, elementfillstatus, results, thresholdfactor, followflowfront)
            [filledelements, newfilledelements, emptyelements] = obj.getFilledElements(elementfillstatus);
            % Default mode: not follow the flowfront
            if nargin < 5
                followflowfront = false;
            end

            highdiffelements = [];
            highdiffsbytimestep = cell(1, size(results, 2));

            for i = 1:size(results,2)
                if followflowfront == true
                    % Treat newly filled elements differently from the rest of
                    % the elements (separate threshold value)
                    oldelements = filledelements{i}(~ismember(filledelements{i}, newfilledelements{i}));
                    threshold = thresholdfactor * mean(nonzeros(results(oldelements, i)));
                    highdiffs = oldelements(results(oldelements, i) > threshold);
                    flowfront_threshold = thresholdfactor * mean(nonzeros(results(newfilledelements{i}, i)));
                    flowfront_highdiffs = newfilledelements{i}(results(newfilledelements{i}, i) > flowfront_threshold);
                    highdiffs_ = [highdiffs; flowfront_highdiffs(~ismember(flowfront_highdiffs, highdiffs))];

                elseif followflowfront == false
                    % Treat all elements the same (one threshold value for
                    % every element) -> default mode
                    threshold = thresholdfactor * mean(nonzeros(results(filledelements{i}, i)));
                    highdiffs = filledelements{i}(results(filledelements{i}, i) > threshold);
                end
                newelements = highdiffs(~ismember(highdiffs, highdiffelements));
                timestepatregister = zeros(size(newelements)) + i;
                newelements = [newelements, timestepatregister];
                highdiffelements = [highdiffelements; newelements];
                highdiffsbytimestep{i} = highdiffs;
            end
            allelementids = 1:size(results, 1);
            normalelements = allelementids(~ismember(allelementids, highdiffelements(:, 1)));
        end
        
        
        function [elementvellen, elementvelangle] = calcElementVelocities(obj, elementnodeids, velx, vely, velz)
            % Using unitvector [1, 0, 0]
            unitvecx = 1;
            unitvecy = 0;
            unitvecz = 0;
            % Length of unitvector
            unitveclen = sqrt(unitvecx.^2 + unitvecy.^2 + unitvecz.^2);
            
            % Velocity components to each node of every element
            elementvelx = obj.getElementValues(elementnodeids, velx);
            elementvely = obj.getElementValues(elementnodeids, vely);
            elementvelz = obj.getElementValues(elementnodeids, velz);
            % Calculate lengths of velocity vectors at each node of every element
            elementvellen = sqrt(elementvelx.^2 + elementvely.^2 + elementvelz.^2);

            % Dotproduct of velocity vector and unitvector at each node every element
            dotproduct = elementvelx.*unitvecx + elementvely.*unitvecy + elementvelz.*unitvecz;
            % Product of the magnitudes of the velocity vectors and unitvectors
            magproduct = elementvellen .* unitveclen;
            % Cosine of the angles of the velocity vectors with the univector
            cosangle = dotproduct ./ magproduct;
            % Angles in radians
            angle_rad = acos(cosangle);
            % Angles in degrees
            elementvelangle = rad2deg(angle_rad);
        end

        function [highdiffs, normaldiffs, highdiffsbytimestep, ideal_threshold, slope, s, frachighdiffs, cutplace_name, fig] = calcHighDiffElementsAutoThreshold(obj, elementfillstatus, elementdiffs, followflowfront, plot_stuff, displayfig, threshstep, threshmax)
            % Setting default values
            if nargin < 8
                threshmax = 5;
            end
        
            if nargin < 7
                threshstep = 0.2;
            end

            if nargin < 6
                displayfig = false;
            end
        
            if nargin < 5
                plot_stuff = false;
            end
        
            if nargin < 4
                followflowfront = false;
            end
        
            % Creating vector with thresholdvalues
            thresholdvalues = transpose(1:threshstep:threshmax);
            highelementdiffs_num = zeros(length(thresholdvalues), 1);
            highelementdiffs_rel = zeros(length(thresholdvalues), 1);
            
            % Calculating number of highdiffs for thresholdvalues
            for i=1:length(thresholdvalues)
                [highelementdiffs, normalelementdiffs] = obj.calcHighDiffElements(elementfillstatus, elementdiffs, thresholdvalues(i), followflowfront);
                highdiffs_num = length(highelementdiffs);
                highdiffs_rel = length(highelementdiffs) / size(elementdiffs, 1);
                 
                % Appending new data
                highelementdiffs_num(i) = highdiffs_num;
                highelementdiffs_rel(i) = highdiffs_rel;
                
                % Calculating cumulative sum and cumulative quadratic sum for
                % variance calculation
                cumsum_1 = cumsum(highelementdiffs_rel(1:i));
                cumsum_2 = cumsum(highelementdiffs_rel(1:i).^2);
                n = transpose(1:length(highelementdiffs_rel(1:i)));
                cumvar = (cumsum_2 - cumsum_1.^2 ./ n) ./ (n-1);
                cumvar(1) = 0;
                
                % Ending loop if zero value or declining cumulative variance is
                % present
                if highdiffs_rel == 0
                    break
                elseif thresholdvalues(i) >= 5 & max(cumvar) ~= cumvar(end)
                    break
                end
            end
            
            % Cutting arrays at loop end
            highelementdiffs_rel = highelementdiffs_rel(1:i);
            highelementdiffs_num = highelementdiffs_num(1:i);
            thresholdvalues = thresholdvalues(1:i);    
            diffs = diff(highelementdiffs_rel);
            [mindiff, mindiffidx] = min(diffs);
            [maxcumvar, maxcumvaridx] = max(cumvar);
            
            cutplace_name = '';
            % Selecting metric to cut the arrays further by
            if diffs(1) - min(diffs) >= 0.3 * (max(diffs) - min(diffs))
                cutat = mindiffidx;
                usemindiff = true;
                usecumvarmax = false;
                cutplace_name = 'deriv min';
            elseif maxcumvar ~= cumvar(end)
                usemindiff = false;
                usecumvarmax = true;
                cutat = maxcumvaridx;
                cutplace_name = 'cumvar max';
            elseif any(highelementdiffs_rel == 0)
                usemindiff = false;
                usecumvarmax = false;
                cutat = find(highdiffelements_rel == 0, 1);
                cutplace_name = 'first 0';
            else
                usemindiff = false;
                usecumvarmax = false;
                cutat = i;
                cutplace_name = 'no cut';
            end
            
            % Cutting the arrays
            highelementdiffs_rel_fit = highelementdiffs_rel(1:cutat);
            thresholdvalues_fit = thresholdvalues(1:cutat);
        
            % Fitting normal distribution to the relative number of highdiffelements
            [m, s] = normfit(highelementdiffs_rel_fit);
            y = normpdf(0:0.01:1, m, s);
            y_true = normpdf(highelementdiffs_rel_fit,m,s);
            [maxval, maxindex] = max(y_true);
            mp = highelementdiffs_rel_fit(maxindex);
            % Fitting linear regression on the thresholdvalues and relative number
            % of highdiffelements
            X = [ones(length(thresholdvalues_fit), 1), thresholdvalues_fit];
            b = X \ highelementdiffs_rel_fit;
            slope = b(2);
            % Calculating ideal threshold
            ideal_threshold = (m - b(1)) / slope;
            
            % Highdiffelements with the ideal threshold
            [highdiffs, normaldiffs, highdiffsbytimestep] = obj.calcHighDiffElements(elementfillstatus, elementdiffs, ideal_threshold, followflowfront);
            numhighdiffs = length(highdiffs);
            frachighdiffs = numhighdiffs / size(elementdiffs, 1);
        
            cutat_disp = thresholdvalues(cutat);
            % Plotting if specified
            if plot_stuff == true
                if displayfig == true
                    fig = figure('Units', 'normalized', 'OuterPosition', [0.1 0.1 0.3 0.8]);
                else
                    fig = figure('Units', 'normalized', 'OuterPosition', [0.1 0.1 0.3 0.8], 'Visible', 'off');
                end
                tiledlayout(5,1);
                ax1 = nexttile;
                scatter(thresholdvalues, highelementdiffs_rel, 'x');
                xline(cutat_disp, 'r');
                title(ax1, 'Rögzített értékek');
                xlabel(ax1, 'Küszöbérték faktor');
                ylabel(ax1, {'Kijelölt elemek'; 'aránya'});
                ax2 = nexttile;
                plot(thresholdvalues(2:end), diffs);
                if usemindiff
                    xline(thresholdvalues(cutat + 1), 'r');
                end
                title(ax2, 'Deriváltak');
                xlabel(ax2, 'Küszöbérték faktor');
                ylabel(ax2, {'Rögzített értékek'; 'deriváltja'});
                ax3 = nexttile;
                plot(thresholdvalues, cumvar);
                if usecumvarmax
                    xline(cutat_disp, 'r');
                end
                title(ax3, 'Kumulált szórásnégyzet');
                xlabel(ax3, 'Küszöbérték faktor');
                ylabel(ax3, {'Kumulált'; 'szórásnégyzet'});
                ax4 = nexttile;
                plot(0:0.01:1, y);
                hold on;
                scatter(highelementdiffs_rel_fit, y_true, 'x')
                title(ax4, 'Az értékekre illesztett normál eloszlás sűrűségfüggvénye');
                text(0.6, max(y)/2, {sprintf('µ = %f', m); sprintf('σ = ±%f', s)});
                xline(m, 'g');
                xline(m - s, 'y');
                xline(m + s, 'y');
                xlabel(ax4, 'Kijelölt elemek aránya');
                ylabel(ax4, {'Sűrűségfüggvény'; 'értéke'});
                ycalc = X*b;
                hold off;
                ax5 = nexttile;
                plot(thresholdvalues_fit, ycalc);
                hold on;
                scatter(thresholdvalues_fit, highelementdiffs_rel_fit, 'x');
                title(ax5, 'Az értékekre illesztett egyenes');
                yline(m, 'g');
                text(ideal_threshold+0.01*(cutat_disp-ideal_threshold), m, ...
                    {sprintf('Meredekség: %f', slope); sprintf('Küszöbfaktor: %f', ideal_threshold)}, 'FontSize',10);
                xline(ideal_threshold, 'g');
                xlabel(ax5, 'Küszöbérték faktor');
                ylabel(ax5, {'Kijelölt elemek'; 'aránya'});
                hold off;
            else
                fig = 0;
            end
        end
        

        % Function for refining the mesh in selected elements by placing multiple
        % nodes
        function [newnodepos, newelementnodeids, affectedelements] = createNewMeshMultiNode(obj, nodepos, elementnodeids, highdiffelements)
            % Create matrix mapping edges to elements
            edges = nchoosek(1:size(elementnodeids, 2), 2);
            edges = transpose(edges);
            edges = reshape(edges, [1, 12]);
            elementedges = elementnodeids(:, edges);
            elementedges = reshape(transpose(elementedges),  2, 6, []);
            elementedges = permute(elementedges, [2, 1, 3]);
            % Choose the edges to refine based on elements chosen for refinement
            edgestorefine = elementedges(:, :, highdiffelements(:, 1));
            edgestorefine = permute(edgestorefine, [2, 1, 3]);
            edgestorefine = transpose(reshape(edgestorefine, 2, []));
            % Remove duplicate rows
            edgestorefine_orig = edgestorefine;
            edgestorefine = unique(sort(edgestorefine, 2), 'rows', 'stable'); % sorting edgestorefine so that edges like 1-2 and 2-1 will be treated as one
            % Get coordinates of the nodes of the edges
            X = nodepos(:, 1);
            Y = nodepos(:, 2);
            Z = nodepos(:, 3);
            refine_X = X(edgestorefine);
            refine_Y = Y(edgestorefine);
            refine_Z = Z(edgestorefine);
            % Calculate the coordinates of the new nodes
            a = 1:length(refine_X);
            a = transpose(a / length(a));
            a = a(randperm(length(a)));
            new_X = mean(refine_X, 2) + 0.1*a.*(refine_X(:, 1) - refine_X(:, 2));
            new_Y = mean(refine_Y, 2) + 0.1*a.*(refine_Y(:, 1) - refine_Y(:, 2));
            new_Z = mean(refine_Z, 2) + 0.1*a.*(refine_Z(:, 1) - refine_Z(:, 2));
            newnodes = [new_X, new_Y, new_Z];
            newnodeids = transpose(size(nodepos, 1)+1:size(nodepos, 1)+size(newnodes, 1));
            % Affected elements -> their edges (sorted) -> newly created nodes ->
            % delaunay-triang.
            elementedges_linear = sort(transpose(reshape(permute(elementedges, [2, 1, 3]), 2, [])), 2);
            affectedelements = unique(ceil(find(ismember(elementedges_linear, edgestorefine, 'rows')) / size(elementedges, 1)));
            
            % Reconnect the vertices and the new points in all affected elements with
            % delaunay-triangulation
            newconnectivity = cell(size(affectedelements, 1), 1);
            for i=1:size(affectedelements, 1)
                % Get hold of the edges and nodes of the current element
                elementedges_current = sort(elementedges(:, :, affectedelements(i)), 2);
                elementnodes_current = unique(reshape(elementedges_current, [], 1));
                elementpoints_current = nodepos(elementnodes_current ,:);
                % Get hold of the respective new nodes to the current element
                newnodeindices_current = find(ismember(edgestorefine, elementedges_current, 'rows'));
                newpoints_current = newnodes(newnodeindices_current, :);
                newnodes_current = newnodeids(newnodeindices_current, :);
                allpoints_current = [elementpoints_current; newpoints_current];
                allnodes_current = [elementnodes_current; newnodes_current];
                DT = delaunayn(allpoints_current);
                %DT_ = delaunay(allpoints_current);
                newconnectivity_current = allnodes_current(DT);
                newconnectivity(i) = {newconnectivity_current};
            end
            
            % Creating the new connectivity matrix
            newconnectivity = cell2mat(newconnectivity);
            newelementnodeids = elementnodeids;
            % Removing refined elements
            newelementnodeids(affectedelements, :) = [];
            % Adding new elements
            newelementnodeids = [newelementnodeids; newconnectivity];
            % Adding new node positions
            newnodepos = [nodepos; newnodes];
        end
        
        
        % Function for calculating the positions of new nodes in each
        % highdiffelement
        function [newnodes, newnodeids] = calcNewNodePos(obj, highdiffelements, values, elementnodeids, nodepos, movecentroid)
            newnodes = [];
            newnodeids = [];
            nodeid = max(elementnodeids, [], "all");
            for i=1:size(highdiffelements, 1)
                nodeid = nodeid + 1;
                elementnodes = elementnodeids(highdiffelements(i, 1), :);
                elementnodecoords = nodepos(elementnodes, :);
                elementvalues = values(highdiffelements(i, 2), elementnodes(:));
                edges = nchoosek(1:length(elementvalues), 2);
                edges_t = transpose(edges);
                pair_values = elementvalues(edges);
                relativeedgediffs = abs(diff(pair_values, 1, 2)) / sum(abs(diff(pair_values, 1, 2)));   
                edgecoords = elementnodecoords(edges_t(:, :), :);
                edgecenters = transpose([mean(reshape(edgecoords(:, 1), 2, []), 1); mean(reshape(edgecoords(:, 2), 2, []), 1); mean(reshape(edgecoords(:, 3), 2, []), 1)]);
                centroid = [mean(elementnodecoords(:, 1)), mean(elementnodecoords(:, 2)), mean(elementnodecoords(:, 3))];
                vectors = edgecenters - centroid;
                sumvector = 0.8 * (sum(vectors .* relativeedgediffs, 1));
                if movecentroid == true
                    newnodepos = centroid + sumvector;
                else
                    newnodepos = centroid;
                end
                newnodes = [newnodes; newnodepos];
                newnodeids = [newnodeids;[elementnodes, nodeid]];
            end
        end
        
        
        % Function for creating the new node connections
        function newnodeconnections = createNewNodeConnections(obj, newnodeids)
            newnodeconnections = [];
            for i=1:size(newnodeids, 1)
                selections = nchoosek(1:length(newnodeids(i, 1:end-1)), size(newnodeids, 2)-2);
                selections = [selections, (zeros(size(selections, 1), 1) + max(selections, [], "all")+1)];
                newconnectivity = reshape(newnodeids(i, selections), size(selections));
                newnodeconnections = [newnodeconnections; newconnectivity];
            end
        end
        
        
        % Function for creating the entire new mesh with node positions and
        % connectivity
        function [newnodepos, newelementnodeids] = createNewMeshOneNode(obj, highdiffelements, values, elementnodeids, nodepos, movecentroid)
            [newnodes, newnodeids] = obj.calcNewNodePos(highdiffelements, values, elementnodeids, nodepos, movecentroid)
            newnodeconnections = obj.createNewNodeConnections(newnodeids)
            % Creating complete new connectivity matrix
            newelementnodeids = elementnodeids
            newelementnodeids(highdiffelements(:, 1), :) = []
            newelementnodeids = [newelementnodeids; newnodeconnections];
            % Appending new node positions to old ones
            newnodepos = [nodepos; newnodes];
        end
    end
end