classdef udmWriter
    % contains only a single funnction for writing udm-s
    
    methods
        % Function for creating the .udm file from the 
        % Function for creating the .udm file from the node positions and
        % connectivity matrix
        function createUDM(newnodepos, newelementnodeids, filename)
            % This block of code generates the txt for the .udm file used for importing
            % the mesh into moldflow.
            
            node_ids = (1:size(newnodepos, 1)).';
            tetra_ids = (1:size(newelementnodeids, 1)).';
            
            % Title block
            title_str = sprintf(['TITL{\n  NAME{"Untitled"}\n  VRSN{"MOLDFLOW"}\n  KEYW{"1713296283"}\n' ...
                                '  DATE{"Apr-16-24"}\n  TIME{"21:38:03"}\n  UDMV{"UDM V3"}\n}']);
            
            % Summary block
            summary_str = sprintf(['SUMR{\n   NOCL{0}\n   NOLY{2}\n   NOTS{1}\n   NOVW{0}\n   NOND{%i}\n' ...
                                   '   NOTX{0}\n   NO1D{0}\n   NOT3{0}\n   NOCV{0}\n   NOSF{0}\n   NORG{0}\n' ...
                                   '   NONB{0}\n   NOSB{0}\n   NOCP{0}\n   NOLC{0}\n   NOT4{%i}\n}'], numel(node_ids), numel(tetra_ids));
            
            % Layer block
            layers_str = sprintf(['LYER{2 " HyperMesh Nodes" 0\n/* NODE */ 4 1 1 0 0\n/* TEXT */ 5 1 1 0 0\n' ...
                                  '/* 1DET */ 6 1 1 0 0\n/* TRI3 */ 7 1 1 0 0\n/* CURV */ 8 1 1 0 0\n/* SURF */ 9 1 1 0 0\n' ...
                                  '/* REGN */ 10 1 1 0 0\n/* NDBC */ 11 1 1 0 0\n/* SUBC */ 12 1 1 0 0\n/* LOCS */ 14 1 1 0 0\n' ...
                                  '/* TET4 */ 15 1 1 0 0\n/* STLR */ 16 1 1 0 0\n}\n' ...
                                  'LYER{4 "vol001" 1\n/* NODE */ 4 1 1 27 0\n/* TEXT */ 5 1 1 27 0\n/* 1DET */ 6 1 1 27 0\n' ...
                                  '/* TRI3 */ 7 1 1 27 0\n/* CURV */ 8 1 1 27 0\n/* SURF */ 9 1 1 27 0\n/* REGN */ 10 1 1 27 0\n' ...
                                  '/* NDBC */ 11 1 1 27 0\n/* SUBC */ 12 1 1 27 0\n/* LOCS */ 14 1 1 27 0\n/* TET4 */ 15 1 1 27 0\n' ...
                                  '/* STLR */ 16 1 1 27 0\n}']);
            
            % Tset block
            tset_str = sprintf(['TSET{40000 1 "Injection location for thermoplastics processes (default)"\n  TCOD{770 ""\n   1.000000e+000}\n' ...
                                '  TCOD{772 ""\n  TCOD{773 ""\n  TCOD{774 ""\n  TCOD{775 ""\n  TCOD{20020 ""}\n  TCOD{20021 ""}\n  TCOD{20023 ""}\n' ...
                                '  TCOD{20032 ""}\n  1.000000e+004 1.000000e+000}\n  TCOD{20040 ""\n  3.001100e+004 1.000000e+000}\n  TCOD{20041 ""\n' ...
                                '  3.001300e+004 1.000000e+000}\n  TCOD{20046 ""\n  3.006000e+004 1.000000e+000}\n  TCOD{20060 ""\n  3.000700e+004 1.000000e+000}\n' ...
                                '  TCOD{20061 ""\n  3.000700e+004 1.000000e+000}\n  TCOD{30503 ""\n  1.000000e+000}\n  TCOD{30504 ""\n  1.000000e+008}\n  TCOD{30505 ""\n' ...
                                '  1.400000e+008}\n  TCOD{30506 ""\n  5.493604e+007}\n  TCOD{30510 ""\n  5.000000e+004 1.000000e+000}\n  TCOD{30520 ""\n  5.001000e+004 1.000000e+000}\n' ...
                                '  TCOD{30530 ""\n  5.002000e+004 1.000000e+000}\n  TCOD{30540 ""\n  5.003000e+004 1.000000e+000}}\n//TSET CARD FOR COOLANT INLETS\n' ...
                                'TSET{40020 1 "Coolant inlet (default) #1"\nTCOD{11102 ""\n  1.666667e-004}\n  TCOD{11103 ""\n  3.000000e+000}\n  TCOD{11104 ""\n' ...
                                '  1.000000e+004}\n  TCOD{11105 ""\n  1.000000e+004}\n  TCOD{11106 ""\n  2.981500e+002}\n  TCOD{20022 ""}}\n TSET{50400 1 "darab  "\n' ...
                                '  TCOD{30107 ""\n        1 }\n  TCOD{925 ""\n        1 }\n  TCOD{11108 ""\n        0 }\n  TCOD{310 ""\n        0 }\n  TCOD{30108 ""\n' ...
	                            '        0 }\n  TCOD{926 ""\n        1 }\n}']);
            
            % Nodal data section
            % Generating nodes_str
            X = newnodepos(:, 1);
            Y = newnodepos(:, 2);
            Z = newnodepos(:, 3);
            
            node_str_base = 'NODE{%d 0 2 1 0.1  %e  %e  %e}';
            node_str_cell = cell(size(newnodepos, 1), 1);
             
            for i = 1:numel(node_ids)
                node_str_cell{i} = sprintf(node_str_base, node_ids(i), X(i), Y(i), Z(i));
            end
            
            nodes_str = strjoin(node_str_cell, newline);
            
            % Tetrahedron data section
            % Generating tetras_str
            node1 = newelementnodeids(:, 1);
            node2 = newelementnodeids(:, 2);
            node3 = newelementnodeids(:, 3);
            node4 = newelementnodeids(:, 4);
            
            tetra_str_base = 'TET4{%d 0 4 27 0 "" 50400 1%s%s%s%s}';
            tetra_str_cell = cell(size(newelementnodeids, 1), 1);
            
            for i = 1:numel(tetra_ids)
                space = ' ';
                node1str = [repelem(space, (8 - numel(num2str(node1(i))))) num2str(node1(i))];
                node2str = [repelem(space, (8 - numel(num2str(node2(i))))) num2str(node2(i))];
                node3str = [repelem(space, (8 - numel(num2str(node3(i))))) num2str(node3(i))];
                node4str = [repelem(space, (8 - numel(num2str(node4(i))))) num2str(node4(i))];
                tetra_str_cell{i} = sprintf(tetra_str_base, tetra_ids(i), node1str, node2str, node3str, node4str);
            end
            
            tetras_str = strjoin(tetra_str_cell, newline);
            
            end_str = sprintf('ENDF');
            
            complete_str = [title_str newline summary_str newline layers_str newline ...
                tset_str newline nodes_str newline tetras_str newline end_str];
            
            fid = fopen(filename, 'w');
            fprintf(fid, complete_str);
            fclose(fid);
        end
    end
end