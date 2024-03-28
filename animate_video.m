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
global cMapInd saveON saveOutputDir exportCAD
cMapInd = 11;% Colourmap index. Check the subfunction "colourGrade" at the end of the file
saveON = false;% Control flag
exportCAD = false;
saveOutputDir='./Fig-EPS-PNG-output';% The output directory
fps = 30;

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

fig = m.fig;
T = 4*fps;

%% Example 1 = Make Mesh for RayleighWave
videoObj = CustomVideoWriter('RayleighWave',saveOutputDir, fps);
videoObj.open(); % Opens the file for writing

for f = 1:T

clc;S=[];% Clear the screen and the function properties
    printProgressBar(f,T,'prefixText','Video Progress', 'id', 'Video','barLength',25);
    fprintf('\n');
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
m.fig = fig;

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

    S.frame = f;
    % Visualization properties
    m.Pplot.limZ(2) = 2;
    S.t = (f*0.5*fps/180)*pi;
    
    visualizationSeqence(m , S)
    
    videoObj.addFrame(m);
end

videoObj.close(); % Closes the video file
%% Example 2 = Make Mesh for Love Wave
videoObj = CustomVideoWriter('LoveWave',saveOutputDir, fps);
videoObj.open(); % Opens the file for writing

for f = 1:T

clc;S=[];% Clear the screen and the function properties
    printProgressBar(f,T,'prefixText','Video Progress', 'id', 'Video','barLength',25);
    fprintf('\n');

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
m.fig = fig;

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

    S.frame = f;
    % Visualization properties
    m.Pplot.limY = [-2 2]*2;
    S.t = (f*0.5*fps/180)*pi;
    
    visualizationSeqence(m , S)

    videoObj.addFrame(m);
end

videoObj.close(); % Closes the video file
%% Example 3 = Make Mesh for P Longitudinal Wave
videoObj = CustomVideoWriter('P_LongitudinalWave',saveOutputDir, fps);
videoObj.open(); % Opens the file for writing

for f = 1:T

clc;S =[];
    printProgressBar(f,T,'prefixText','Video Progress', 'id', 'Video','barLength',25);
    fprintf('\n');

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
m.fig = fig;

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


    S.frame = f;
    % Visualization properties
    m.Pplot.limZ(2) = 2;
    S.t = (f*0.5*fps/180)*pi;
    
    visualizationSeqence(m , S)

    videoObj.addFrame(m);
end

videoObj.close(); % Closes the video file

%% Example 4 = Make Mesh for S Transverse Wave
videoObj = CustomVideoWriter('S_TransverseWave',saveOutputDir, fps);
videoObj.open(); % Opens the file for writing

for f = 1:T


clc;S =[];
    printProgressBar(f,T,'prefixText','Video Progress', 'id', 'Video','barLength',25);
    fprintf('\n');

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
m.fig = fig;    

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


    S.frame = f;
    % Visualization properties
    m.Pplot.limZ(2) = 2;
    S.t = (f*0.5*fps/180)*pi;
    
    visualizationSeqence(m , S)

    videoObj.addFrame(m);
end

videoObj.close(); % Closes the video file