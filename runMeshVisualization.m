% 
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

% function runMeshVisualization(cMapInd)

%% Initialize Mesh Generator

% Clean the workspace
clc;close all;

% converFig2EPS_PNG_PDF_submodule
addpath("./converFig2EPS_PNG_PDF_submodule")
% Add the displacement/deformation functions
addpath("./Displacement&Deformation_Functions");

% Properties
cMapInd = 11;% Colourmap index. Check the subfunction "colourGrade" at the end of the file
saveON = false;% Control flag
saveOutputDir='./Fig-EPS-PNG-output';% The output directory

% Initialize the mesh generator
m = MeshGenerator();
% Setup the mesh
m.setupMesh();
% Plot the cube
m.plotElements(m.translateElements(m.cube,[-1 -1 -1]*0.5*0));

% Get the xyz mesh dimentions and properties
xyz = m.xyz;

%% Example 1 = Make Mesh for RayleighWave
clc;S=[];% Clear the screen and the function properties

    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention index. Very easy way to do dimention bisection
    xyz.N = [50 3 20]';% Number of mesh elements
    xyz.gapP = [20 20 20]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-10 10];% Limits
    xyz{'y','Lim'} = [-1.5 1.5];% Limits
    xyz{'z','Lim'} = [-3 0];% Limits
    
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor
    m.cube.FaceAlpha = 100/100;

% Make mesh
m.setupMesh(xyz);

% Plot Mesh
m.plotMeshFast(m.translateElements(m.Elements,[0 0 0]*m.Pplot.scaleFactor));

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function 
    S.i = 1;% Proerpty set index 
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


% Visualization properties
m.Pplot.limZ(2) = 2;
% m.Pplot.axisEqual = 1;
% m.Pplot.View = 3;
m.cube.colourF = [204 204 0]/255;
m.cube.colourE = 70/100*[1 1 1];%[140 140 0]/255;
m.cube.colourE = [100 100 100]/255;


% Get elements set "EE" which is the set we will manipulate
% This is actually the mesh as an element set
EE = m.Elements;

% Transformation: Translate the mesh
% EE = m.translateElements(EE,[0 0 0]*m.Pplot.scaleFactor);

% Transformation: Deformation or displacement
switch (1)
    case 1% Deformation plotting        
        EE = m.transformNodesR(EE,@(R) RayleighWave_Displacements(R,S));
    case 2% Translation plotting
        EE = m.transformElements(EE,@(E) displaceCenter_Uxyz(E,S));
end    

% Update colours
EE = m.transformElements(EE,@(E) colourUpdate(E,m.cube));
% Calculate the colour gradients
R0 = m.getMatrix(m.Elements,'nodesR');
R  = m.getMatrix(EE,'nodesR');
distM = max(sqrt(sum((R - R0).^2,2)));
EE = m.transformElements(EE,@(E) colourGrade(E,(m.Elements),distM,cMapInd));%struct2table
    
% Ploting element by element or all using the matrix. 
% The difference is that the fast method doesn't allow for custom colour
% grading while the slow method being element by element allows for percise
% coloring
if cMapInd<0; m.plotMeshFast(EE); else
m.plotElements(EE); end

% Post polot adjustments
axis off;

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


fileNames{S.i} = 'RayleighWave';
if saveON
    fileNames{S.i} = func2str(S.Uxyz);
    fileNames{S.i} = [fileNames{S.i},'_C_',num2str(cMapInd)];
    savefig(m.fig,fullfile(saveOutputDir,fileNames{S.i}));
end

%% Example 2 = Make Mesh for Love Wave
clc;S=[];% Clear the screen and the function properties
    
    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention
    xyz.N = [100 3 10]';% Number of mesh points
    xyz.gapP = [0 20 0]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-8 8];% Limits
    xyz{'y','Lim'} = [-2 2];% Limits
    xyz{'z','Lim'} = [-4 0];% Limits
    
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor
    m.cube.FaceAlpha = 100/100;

% Make mesh
m.setupMesh(xyz);

% Plot Mesh
% m.plotMeshFast(m.translateElements(m.Elements,[0 0 0]*m.Pplot.scaleFactor));

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function 
    S.i = 2;
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

% View properties
m.Pplot.limY = [-2 2]*2;
% m.Pplot.axisEqual = 1;
% m.Pplot.View = 3;
m.cube.colourF = [204 204 0]/255;
m.cube.colourE = 70/100*[1 1 1];%[140 140 0]/255;
m.cube.colourE = [100 100 100]/255;

%  Get elements 
EE = m.Elements;
% Translate
EE = m.translateElements(EE,[0 0 0]*m.Pplot.scaleFactor);
% Make Displacements
switch (1)
    case 1% Deformation plotting        
        EE = m.transformNodesR(EE,@(R) LoveWave_Displacements(R,S));
    case 2% Translation plotting
        EE = m.transformElements(EE,@(E) displaceCenter_Uxyz(E,S));
end    
% Update colours
EE = m.transformElements(EE,@(E) colourUpdate(E,m.cube));
% Calculate the colour gradients
R0 = m.getMatrix(m.Elements,'nodesR');
R  = m.getMatrix(EE,'nodesR');
distM = max(sqrt(sum((R - R0).^2,2)));
EE = m.transformElements(EE,@(E) colourGrade(E,(m.Elements),distM,cMapInd));%struct2table
    
% Plot 
if cMapInd<0; m.plotMeshFast(EE); else
m.plotElements(EE); end

% Post polot adjustments
axis off
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

fileNames{S.i} = 'LoveWave';
if saveON
    fileNames{S.i} = func2str(S.Uxyz);
    fileNames{S.i} = [fileNames{S.i},'_C_',num2str(cMapInd)];
    savefig(m.fig,fullfile(saveOutputDir,fileNames{S.i}));
end

%% Example 3 = Make Mesh for P Longitudinal Wave
clc;S =[];

    % Mesh Properties 
    m.dimRemoveID = 0;% Remove a dimention
    xyz.N = [150 3 10]';% Number of mesh points
    xyz.gapP = [0 20 0]'./100;% Gap percetnage 
    xyz{'x','Lim'} = [-8 8];% Limits
    xyz{'y','Lim'} = [-1.5 1.5];% Limits
    xyz{'z','Lim'} = [-3 0];% Limits
    
    % Plotting properties
    m.Pplot.scaleFactor = 15/10;% Axes scale factor
    m.cube.FaceAlpha = 100/100;

% Make mesh
m.setupMesh(xyz);

% Plot Mesh
% m.plotMeshFast(m.translateElements(m.Elements,[0 0 0]*m.Pplot.scaleFactor));

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function
    S.i = 3;
    S.Uxyz = @P_LongitudinalWave_Displacements;
    S.t = (0*10/180)*pi;% Time
    % Material Properties for the deformation function: P Wave
    S.A = 5/10;
    S.k = 10/10;
    S.omega = 1;

% View properties
m.Pplot.limZ(2) = 2;
% m.Pplot.axisEqual = 1;
% m.Pplot.View = 3;
m.cube.colourF = [204 204 0]/255;
m.cube.colourE = 70/100*[1 1 1];%[140 140 0]/255;
m.cube.colourE = [100 100 100]/255;


%  Get elements 
EE = m.Elements;
% Translate
% EE = m.translateElements(EE,[0 0 0]*m.Pplot.scaleFactor);
% Make Displacements
switch (1)
    case 1% Deformation plotting        
        EE = m.transformNodesR(EE,@(R) P_LongitudinalWave_Displacements(R,S));
    case 2% Translation plotting
        EE = m.transformElements(EE,@(E) displaceCenter_Uxyz(E,S));
end    
% Update colours
EE = m.transformElements(EE,@(E) colourUpdate(E,m.cube));
% Calculate the colour gradients
R0 = m.getMatrix(m.Elements,'nodesR');
R  = m.getMatrix(EE,'nodesR');
distM = max(sqrt(sum((R - R0).^2,2)));
EE = m.transformElements(EE,@(E) colourGrade(E,(m.Elements),distM,cMapInd));%struct2table
    
% Plot 
if cMapInd<0; m.plotMeshFast(EE); else
m.plotElements(EE); end

% Post polot adjustments
axis off
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

fileNames{S.i} = 'P_LongitudinalWave';
if saveON
    fileNames{S.i} = func2str(S.Uxyz);
    fileNames{S.i} = [fileNames{S.i},'_C_',num2str(cMapInd)];
    savefig(m.fig,fullfile(saveOutputDir,fileNames{S.i}));
end

%% Example 4 = Make Mesh for S Transverse Wave
clc;S =[];
% Mesh Properties 
m.dimRemoveID = 0;% Remove a dimention
xyz.N = [150 3 20]';% Number of mesh points
xyz.gapP = [0 20 20]'./100;% Gap percetnage 
xyz{'x','Lim'} = [-10 10];% Limits
xyz{'y','Lim'} = [-1.5 1.5];% Limits
xyz{'z','Lim'} = [-3 0];% Limits

% Plotting properties
m.Pplot.scaleFactor = 15/10;% Axes scale factor
m.cube.FaceAlpha = 100/100;

% Make mesh
m.setupMesh(xyz);

% Plot Mesh
% m.plotMeshFast(m.translateElements(m.Elements,[0 0 0]*m.Pplot.scaleFactor));

    % Diplay Wave settings determinted by the variables defined in the "S.Uxyz" referenced function
    S.i = 4;
    S.Uxyz = @S_TransverseWave_Displacements;
    S.t = (0*10/180)*pi;% Time
    % Material Properties for Rayleigh Wave
    S.k = 10/10;
    S.A = 10/10;
    S.omega = 1;


% View properties
m.Pplot.limZ(2) = 2;
% m.Pplot.axisEqual = 1;
% m.Pplot.View = 3;
m.cube.colourF = [204 204 0]/255;
m.cube.colourE = 70/100*[1 1 1];%[140 140 0]/255;
m.cube.colourE = [100 100 100]/255;


%  Get elements 
EE = m.Elements;
% Translate
% EE = m.translateElements(EE,[0 0 0]*m.Pplot.scaleFactor);
% Make Displacements
switch (1)
    case 1% Deformation plotting        
        EE = m.transformNodesR(EE,@(R) S_TransverseWave_Displacements(R,S));
    case 2% Translation plotting
        EE = m.transformElements(EE,@(E) displaceCenter_Uxyz(E,S));
end    
% Update colours
EE = m.transformElements(EE,@(E) colourUpdate(E,m.cube));
% Calculate the colour gradients
R0 = m.getMatrix(m.Elements,'nodesR');
R  = m.getMatrix(EE,'nodesR');
distM = max(sqrt(sum((R - R0).^2,2)));
EE = m.transformElements(EE,@(E) colourGrade(E,(m.Elements),distM,cMapInd));%struct2table
    
% Plot 
if cMapInd<0; m.plotMeshFast(EE); else
m.plotElements(EE); end

% Post polot adjustments
axis off
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

fileNames{S.i} = 'S_TransverseWave';
if saveON
    fileNames{S.i} = func2str(S.Uxyz);
    fileNames{S.i} = [fileNames{S.i},'_C_',num2str(cMapInd)];
    savefig(m.fig,fullfile(saveOutputDir,fileNames{S.i}));
end

%% Convert Fig to EPS
clc;close all;
% Convert FIG to: PNG, EPS, PDF
convertFig2Eps(saveOutputDir,'png','expand','expGraph','eval: ax=fig.Children;  objOutH=ax; ResolutionDPI=450;');
convertFig2Eps(saveOutputDir,'eps','expand','expGraph','eval: ax=fig.Children;  objOutH=ax; ResolutionDPI=450;');
convertFig2Eps(saveOutputDir,'','expand','expGraph','eval: ax=fig.Children;  objOutH=ax; ResolutionDPI=450; outputExt=''.pdf'';');
close all; 

% end

%% Sub functions 
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
%     E.FaceAlpha=1;
end

% Get the colours from the cube and pass it to the element
function E = colourUpdate(E,cube)
% function E = colourUpdate(E,cube)
    E.colourF = cube.colourF;
    E.colourE = cube.colourE;
    E.FaceAlpha = cube.FaceAlpha;
end