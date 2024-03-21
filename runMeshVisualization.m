% 
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

%% Initialize Mesh Generator

% Clean the workspace
clear; clc;close all; 

% converFig2EPS_PNG_PDF_submodule
addpath("./converFig2EPS_PNG_PDF_submodule")
% Add the displacement/deformation functions
addpath("./Displacement&Deformation_Functions");

% Properties
global cMapInd saveON saveOutputDir
cMapInd = 11;% Colourmap index. Check the subfunction "colourGrade" at the end of the file
saveON = false;% Control flag
saveOutputDir='./Fig-EPS-PNG-output';% The output directory

% Initialize the mesh generator
m = MeshGenerator();
% Setup the mesh
m.setupMesh();
% Plot the cube
m.plotElements(m.translate(m.cube,[-1 -1 -1]*0.5*0));

% Get the xyz mesh dimentions and properties
xyz = m.xyz;

% Plotting properties
m.Pplot.scaleFactor = 15/10;% Axes scale factor

%% Example 1 = Make Mesh for RayleighWave
clc;S=[];% Clear the screen and the function properties

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function 
    S.i = 1;% Proerpty set index 
    S.fileName = 'RayleighWave';
    S.Uxyz = @RayleighWave_Displacements;% Deformation function handle/reference
    % Function input time
    S.t = (0*10/180)*pi;
    % Material Properties for the deformation function: Rayleigh Wave
    S.c = 98/100;
    S.cl = 85/100;
    S.ct = 74/100;
    S.A1 = 34/10;
    S.A2 = 34/10;
    S.omega = 1;

% Initialize the mesh generator
m = MeshGenerator();

    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention index. Very easy way to do dimention bisection
    xyz.N = [50 3 20]';% Number of mesh elements
    xyz.gapP = [20 20 20]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-10 10];% Limits
    xyz{'y','Lim'} = [-1.5 1.5];% Limits
    xyz{'z','Lim'} = [-3 0];% Limits   
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor

% Make mesh
m.setupMesh(xyz);


% Visualization properties
m.Pplot.limZ(2) = 2;

visualizationSeqence(m , S)

%% Example 2 = Make Mesh for Love Wave
clc;S=[];% Clear the screen and the function properties

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function 
    S.i = 2;
    S.fileName = 'LoveWave';
    S.Uxyz = @LoveWave_Displacements;
    S.t = (0*10/180)*pi;% Time
    % Material Properties for the deformation function: Love Wave
    S.c = 87/100;   
    S.c1 = 74/100;
    S.c2 = 95/100;
    S.A1 = 1/10;
    S.A2 = 20/10;
    S.h  = -1;
    S.omega = 10/10;

% Initialize the mesh generator
m = MeshGenerator();

    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention
    xyz.N = [100 3 10]';% Number of mesh points
    xyz.gapP = [0 20 0]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-8 8];% Limits
    xyz{'y','Lim'} = [-2 2];% Limits
    xyz{'z','Lim'} = [-4 0];% Limits
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor
    
% Make mesh
m.setupMesh(xyz);

% View properties
m.Pplot.limY = [-2 2]*2;

visualizationSeqence(m , S)

%% Example 3 = Make Mesh for P Longitudinal Wave
clc;S =[];

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function
    S.i = 3;
    S.fileName = 'P_LongitudinalWave';
    S.Uxyz = @P_LongitudinalWave_Displacements;
    S.t = (0*10/180)*pi;% Time
    % Material Properties for the deformation function: P Wave
    S.A = 5/10;
    S.k = 10/10;
    S.omega = 1;

% Initialize the mesh generator
m = MeshGenerator();

    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention
    xyz.N = [150 3 10]';% Number of mesh points
    xyz.gapP = [0 20 0]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-8 8];% Limits
    xyz{'y','Lim'} = [-1.5 1.5];% Limits
    xyz{'z','Lim'} = [-3 0];% Limits
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor

% Make mesh
m.setupMesh(xyz);


% View properties
m.Pplot.limZ(2) = 2;

visualizationSeqence(m , S)

%% Example 4 = Make Mesh for S Transverse Wave
clc;S =[];

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function
    S.i = 4;
    S.fileName = 'S_TransverseWave';
    S.Uxyz = @S_TransverseWave_Displacements;
    S.t = (0*10/180)*pi;% Time
    % Material Properties for Rayleigh Wave
    S.k = 10/10;
    S.A = 10/10;
    S.omega = 1;

% Initialize the mesh generator
m = MeshGenerator();
    
    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention
    xyz.N = [150 3 20]';% Number of mesh points
    xyz.gapP = [0 20 20]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-10 10];% Limits
    xyz{'y','Lim'} = [-1.5 1.5];% Limits
    xyz{'z','Lim'} = [-3 0];% Limits
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor
    
% Make mesh
m.setupMesh(xyz);

% View properties
m.Pplot.limZ(2) = 2;

visualizationSeqence(m , S)

%% Convert Fig to EPS
if saveON
    clc;close all;
    % Convert FIG to: PNG, EPS, PDF
    convertFig2Eps(saveOutputDir,'png','expand','expGraph','eval: ax=fig.Children;  objOutH=ax; ResolutionDPI=450;');
    convertFig2Eps(saveOutputDir,'eps','expand','expGraph','eval: ax=fig.Children;  objOutH=ax; ResolutionDPI=450;');
    convertFig2Eps(saveOutputDir,'','expand','expGraph','eval: ax=fig.Children;  objOutH=ax; ResolutionDPI=450; outputExt=''.pdf'';');
    close all; 
end

%% Sub functions 
function visualizationSeqence(m , S)
    global cMapInd saveON saveOutputDir

    % Mesh properties
    % m.Pplot.axisEqual = true;
    % m.Pplot.View = 3;
    m.cube.colourF = [204 204 0]/255;
    m.cube.colourE = 70/100*[1 1 1];%[140 140 0]/255;
    m.cube.colourE = [100 100 100]/255;
    m.cube.faceAlpha = 100/100;

    % Plot Mesh
    m.plotMeshFast(m.translate(m.elements,[0 0 0]*m.Pplot.scaleFactor));

    % Get elements set "EE" which is the set we will manipulate
    % This is actually the mesh as an element set
    EE = m.elements;
    EE = m.transformElements(EE,@(E) E.updateColourPropertiesFrom(m.cube));
    R0 = m.getMatrix(EE,'nodesR'); % Get initial state
    EE0 = m.copyElements();
    
    % Transformation: Translate the mesh
    % EE = m.translate(EE,[0 0 0]*m.Pplot.scaleFactor);
    
    % Transformation: Deformation or displacement
    switch (1)
        case 1% Deformation plotting        
            EE = m.transformNodesR(EE,@(R) S.Uxyz(R,S));
        case 2% Translation plotting
            EE = m.transformElements(EE,@(E) displaceCenter_Uxyz(E,S));
    end
    % Get matrix state post process
    R  = m.getMatrix(EE,'nodesR');
    
    % Calculate the colour gradients    
    distM = max(sqrt(sum((R - R0).^2,2)));
    EE = m.transformElements(EE,@(E) colourGrade(E,(EE0),distM,cMapInd));%struct2table
        
    % Ploting element by element or all using the matrix. 
    % The difference is that the fast method doesn't allow for custom colour
    % grading while the slow method being element by element allows for percise
    % coloring
    if cMapInd<0; m.plotMeshFast(EE); else
    m.plotElements(EE); end
    
    % Post polot adjustments
    axis off;
    % m.Pplot.axP.LineWidth = 0.0001;
    % set(get(m.fig.Children.Children,'children'), 'edgecolor', [0 0 0])
    
    % Rotation and view -- if needed
    % rotate(m.fig.Children.Children,[1 0 0],90,mean(xyz.Lim,2));
    % rotate(m.Pplot.axP,[1 0 0],90,mean(xyz.Lim,2))
    % view(157,-66)
    % camproj('orthographic');
    % camproj('perspective');
    % set(gca, 'CameraPosition', [0,0,0])
    % set(gca, 'CameraTarget', [0,0,1])
    % set(gca, 'CameraUpVector', [0,1,1])
    
    if saveON
        S.fileName = func2str(S.Uxyz);
        S.fileName = [S.fileName,'_C_',num2str(cMapInd)];
        savefig(m.fig,fullfile(saveOutputDir,S.fileName));
    end
    disp("Done!")
    m.exportToFile(saveOutputDir+"/"+S.fileName,"STL")
    m.exportToFile(saveOutputDir+"/"+S.fileName,"OBJ")

end

function E = displaceCenter_Uxyz(E,S)
    rc = E.center;
    T = S.Uxyz(rc,S)-rc;
    E.nodesR = E.nodesR + (10/10)*repmat(T,size(E.nodesR,1),1);
end
            
% This function is used for colour grading
function E = colourGrade(E,Es_original,distM,cmapInd)
% Es_original_TBL
%     E0 = table2struct(Es_original_TBL(Es_original_TBL.id == E.id,:));
    E0 = Es_original(E.id);
    
    R0 = E0.('nodesR');
    R  = E.('nodesR');    
    p = mean(sqrt(sum((R - R0).^2,2)))./distM;
    
   cmaps = {'parula','turbo','hsv','hot','cool','spring','summer',...
            'autumn','winter','gray','bone', 'copper','pink', 'lines','jet',...
            'colorcube','prism','flag'};
    if 0 < cmapInd && cmapInd <= numel(cmaps)
        E.colourF = vals2colormap(p,cmaps{cmapInd},[0 1]);    
    else
        E.colourF = [238, 104, 104]/255*p + [204 204 0]/255*(1-p);
    end
    E.faceAlpha=p;
    E.edgeAlpha=p;
    E.colourE = "none"; %(1- E.colourF) *(1-E.faceAlpha); 
%     E.nodesR = gpuArray(E.nodesR);
%     E.FaceConnectNodes = gpuArray(E.FaceConnectNodes);
end
