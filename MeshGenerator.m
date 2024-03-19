%
% Class-name:   MeshGenerator 
% Description: This script generates 2D/3D mesh and visualizes it.
% 
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

%% Mesh generator
classdef MeshGenerator < handle
    
    % Internal properties
    properties (SetAccess=private)
        cubeTemplate;% This is a cube prototype
        elementAlignemntTypes;
        Elements;% Store the elements
        xyz;% The dimension and setup of the mesh
    end

    % Public properties
    properties (SetAccess=public)
        cube;% Working cube modified  
    end
    properties (SetAccess=public)
        dimRemoveID;
        Pplot; 
        fig;
    end
    
    methods
        % Constructor
        function [mesh] = MeshGenerator()
            % Settings
            mesh.dimRemoveID = 0;% 0=>OFF 1=>x 2=>y 3=>z
            % Labels
            xyz(1).Label='x';
            xyz(2).Label='y';
            xyz(3).Label='z';
            for i = 1:numel(xyz)
                % Constant IDs
                xyz(i).ID=i;
                % Mesh properties: Default Limits
                xyz(i).Lim = [0 1]; 
                % Mesh properties: Default Number of Elements
                xyz(i).N = 1;
                % Element Gap factor between elements
                xyz(i).gapP= 0.1;% Input 
                % Element side length
                xyz(i).ElementL = mesh.calcElSpacing(xyz(i))*(1-xyz(i).gapP); 
                % Element side length gap
                xyz(i).gap =  mesh.calcElSpacing(xyz(i))*(xyz(i).gapP);
                % Mesh grid coordinates
                xyz(i).Coords = {};
                % Mesh grid coordinates
                xyz(i).mesh =[];
            end

            % Element Alignment
            mesh.elementAlignemntTypes = ["Symetric";"Center"];
            % Default Plot view direction
            mesh.Pplot.View = 3;
            mesh.Pplot.gridON = 1;
            mesh.Pplot.scaleFactor = 1;
            mesh.Pplot.axisEqual = 1;
            % Cube Template Coordinates
            mesh.cube.id = 0;
            mesh.cube.nodesR=[0 0 0; 0 0 1; 1 0 1; 1 0 0; 0 1 0; 0 1 1; 1 1 1; 1 1 0; ];
            mesh.cube.FaceConnectNodes = [1 4 3 2; 5 8 7 6; 1 2 6 5; 3 4 8 7; 2 3 7 6; 1 5 8 4];
            mesh.cube.center = mean(mesh.cube.nodesR);
            % Properties
            mesh.cube.colourF = [100 100 100]/100;%RGB percentage or use 'none'
            mesh.cube.colourE = [0 0 0]/100;
            mesh.cube.FaceAlpha = 100/100;% Alpha percentage
            
            % Assign Mesh
            mesh.xyz = struct2table(xyz,'RowNames',{xyz.Label});
            mesh.cubeTemplate = mesh.cube;% Record the original type
        end
        
        % Setup mesh
        function setupMesh(mesh,xyz_)% mesh obj
            % Update the mesh parameters 
            if nargin == 2 && istable(xyz_)
                mesh.xyz = xyz_;
%                 disp(mesh.xyz);
                disp('Table updated!');
            else
                disp('Table not updated!');
            end
            
            % Setup Properties and Update colours
            mesh.cubeTemplate.colourF    = mesh.cube.colourF;
            mesh.cubeTemplate.colourE    = mesh.cube.colourE;
            mesh.cubeTemplate.FaceAlpha = mesh.cube.FaceAlpha;
            % Update the colours of the template cube
            mesh.cube = mesh.cubeTemplate;
            
            % Default viewing Limits 
            mesh.Pplot.limX = mesh.xyz{'x','Lim'}*mesh.Pplot.scaleFactor;
            mesh.Pplot.limY = mesh.xyz{'y','Lim'}*mesh.Pplot.scaleFactor;
            mesh.Pplot.limZ = mesh.xyz{'z','Lim'}*mesh.Pplot.scaleFactor;                       
            
            % Make coordinates
            mesh.xyz('x',:) = mesh.setupCoordinate(mesh.xyz('x',:));
            mesh.xyz('y',:) = mesh.setupCoordinate(mesh.xyz('y',:));
            mesh.xyz('z',:) = mesh.setupCoordinate(mesh.xyz('z',:));

            % Remove the dimension component 
            if mesh.dimRemoveID > 0 && mesh.dimRemoveID <= 3
                % Remove the cube nodes that correspond to the removed dimension   
                ind = not(mesh.cube.nodesR(:,mesh.dimRemoveID) > 0);
                mesh.cube.nodesR = mesh.cube.nodesR(ind,:);

                % Get the right square face
                % For removing For X
                if mesh.dimRemoveID == 1                    
                    mesh.cube.FaceConnectNodes = [1 2 4 3];
                end
                % For removing For Y
                if mesh.dimRemoveID == 2                    
                    mesh.cube.FaceConnectNodes = [1 2 3 4];
                end
                % For removing For Z
                if mesh.dimRemoveID == 3
                    mesh.cube.FaceConnectNodes = [1 2 4 3];
                end
            end
            
            % Scale the cube
            mesh.cube = mesh.scaleElements(mesh.cube,...
                [mesh.calcElSpacing(mesh.xyz('x',:))  ...
                 mesh.calcElSpacing(mesh.xyz('y',:)) ...
                 mesh.calcElSpacing(mesh.xyz('z',:))].*(1-mesh.xyz.gapP'));

            % Cube Test preview
%             mesh.PplotElements(mesh.translateElements(mesh.cube,[-1 -1 -1]*0.5*0));
            
            [X,Y,Z] = meshgrid(mesh.xyz{"x","Coords"}{:},mesh.xyz{"y","Coords"}{:},mesh.xyz{"z","Coords"}{:});
            % Assign Mesh
            mesh.xyz{"x","mesh"} = {X};
            mesh.xyz{"y","mesh"} = {Y};
            mesh.xyz{"z","mesh"} = {Z};
            
            Elements_ = [];
            % Make mesh
            for i = 1:numel(X)
                % Make elements
                El = mesh.cube;
                El.id = i;
                % Translate the 
                El = mesh.translateElements(El,[X(i),Y(i),Z(i)]);
                % For fast plotting
                El.FaceConnectNodes_Fast = El.FaceConnectNodes + (i-1)*(8*(mesh.dimRemoveID==0)+4*(mesh.dimRemoveID~=0)); 
%                 El.colourF_Fast     = repmat(El.colourF,   size(El.FaceConnectNodes_Fast,1),1);
%                 El.FaceAlpha_Fast  = repmat(El.FaceAlpha,size(El.FaceConnectNodes_Fast,1),1);
                % Gather elements together
                Elements_ = [Elements_;El];
            end
            mesh.Elements = Elements_;
            % Plot Preview
%             mesh.PplotElements(mesh.translateElements(Elements_,[0 0 0]*0.5*0));


        end
                
        % Setup Coordinate
        function [c] = setupCoordinate(mesh,c)%s - general coordinates, P - properties

            % Setup Mesh properties: Element length
            c.ElementL = mesh.calcElSpacing(c)*(1-c.gapP);
            % Setup Mesh properties: Element gap
            c.gap      = mesh.calcElSpacing(c)*(c.gapP);

            % Remove a dimension for the grid
            c.Lim = c.Lim*not(mesh.dimRemoveID==c.ID);

            % Make a grid origins i.e. the [0 0 0] vertex of the Element cube polygon
            Coords = linspace(c.Lim(1),c.Lim(2),(c.N+1));
            Coords(end) =[]; 
            c.Coords = {Coords};
        end
        
        function size_element = calcElSpacing(~,c)
        % function size_element = sizeElement(~,c)
            size_element = (c.Lim(2)-c.Lim(1))./c.N;
        end
            
        % Subfunction for plotting element
        function plotElements(mesh,elements)
            % Get handles and do plotting
            try 
                if not(isvalid(mesh.fig))
                    error();
                end
            catch
                mesh.fig = figure();
                mesh.fig.WindowStyle = 'docked';
            end
            clf(mesh.fig);
            view(mesh.Pplot.View);
            ax = mesh.fig.Children(1);
            % Get axes and setup plot limits
%             ax = mesh.ax;
            if mesh.Pplot.axisEqual
                axis(ax,'equal');
            end
            ax.XLim = mesh.Pplot.limX;
            ax.YLim = mesh.Pplot.limY;
            ax.ZLim = mesh.Pplot.limZ;
            ax.XLabel.String = 'x';
            ax.YLabel.String = 'y';
            ax.ZLabel.String = 'z';
            if mesh.Pplot.gridON
                grid(ax,'on')
            end
            mesh.Pplot.axP = [];
            % Plot elements one by one
            for i = 1:numel(elements)
                E = elements(i);
                mesh.Pplot.axP(i) = patch(ax,'Faces', E.FaceConnectNodes, 'Vertices', E.nodesR, ...
                'Facecolor', E.colourF,'FaceAlpha',E.FaceAlpha,'EdgeColor',E.colourE); 
        %         drawnow
            end
            drawnow;
        end

        % Subfunction for plotting element
        function plotMeshFast(mesh,elements)
            % Get handles and do plotting
            try 
                if not(isvalid(mesh.fig))
                    error;
                end
            catch
                mesh.fig = figure();
                mesh.fig.WindowStyle = 'docked';
            end
            clf(mesh.fig);
            view(mesh.Pplot.View);
            ax = mesh.fig.Children(1);
            % Get axes and setup plot limits
%             ax = mesh.ax;
            if mesh.Pplot.axisEqual
                axis(ax,'equal');
            end
            ax.XLim = mesh.Pplot.limX;
            ax.YLim = mesh.Pplot.limY;
            ax.ZLim = mesh.Pplot.limZ;
            ax.XLabel.String = 'x';
            ax.YLabel.String = 'y';
            ax.ZLabel.String = 'z';
            if mesh.Pplot.gridON
                grid(ax,'on')
            end
            mesh.Pplot.axP = [];
            
            % Get properties ax matrix
            FaceConnectNodes_Fast = mesh.getMatrix(elements,"FaceConnectNodes_Fast");            
            nodesR = mesh.getMatrix(elements,"nodesR");
%             colour =  mesh.getMatrix(elements,"Colour_Fast");
%             FaceAlpha = mesh.getMatrix(elements,"FaceAlpha_Fast");
            
            mesh.Pplot.axP = ...
            patch(ax,'Faces', FaceConnectNodes_Fast, 'Vertices', nodesR,...
            'FaceColor',mesh.cube.colourF,'FaceAlpha',mesh.cube.FaceAlpha, 'EdgeColor', mesh.cube.colourE); 
            drawnow;
        end
        
        % Subfunction for translating element
        function [elements] = translateElements(~,elements,translationVector)
            for i = 1:numel(elements)
                % Translate vector        
                elements(i).nodesR = elements(i).nodesR + repmat(translationVector,size(elements(i).nodesR,1),1);
                % Calculate the center
                elements(i).center   = mean(elements(i).nodesR);
            end
        end

        % Subfunction for scaling element
        function [Elements] = scaleElements(~,Elements,scale)
            for i = 1:numel(Elements)
                R = Elements(i).nodesR;
                c = Elements(i).center;
                % Scale vector        
                Elements(i).nodesR = R.*scale;

                % Calculate the center
                Elements(i).center   = mean(Elements(i).nodesR);
            end
        end
        
        % Transformation functions
        function [elements] = transformElements(~,elements,transformationF)
        % function [elements] = transformElements(~,elements,transformationF)   
            for i = 1:numel(elements)
                % Translate vector        
                elements(i) = transformationF(elements(i));                
            end
        end
        
        function [elements] = transformNodesR(~,elements,transformationF)
        % function [elements] = transformNodesR(~,elements,transformationF)   
            for i = 1:numel(elements)
                % Transform Nodes
                elements(i).nodesR = transformationF(elements(i).nodesR);                
            end
        end
        
        
        %%%%% Add function for element rotation
        
        % Get method for the max XYZ
        function xyz = get.xyz(mesh)
             xyz = mesh.xyz;
        end
        
        % Convert elements to a matrix of coordinates
        function M = getMatrix(~,elements,property)
        % function M = getMatrix(~,elements,property)
        % M = cell2mat({elements.(property)}');
            M = cell2mat({elements.(property)}');
        end
    end
end


