#
# Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
#
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mesh_generator import MeshGenerator
from colour_grading import vals2colormap, displace_center_uxyz, colour_grade, colour_update
from displacement_and_deformation_functions import (
    love_wave_displacements,
    p_longitudinal_wave_displacements,
    rayleigh_wave_displacements,
    s_transverse_wave_displacements,
)


# function runMeshVisualization(cMapInd)

## --- Initialize Mesh Generator --- 

# Clear all figures
plt.close('all')

# converFig2EPS_PNG_PDF_submodule
sys.path.append("./converFig2EPS_PNG_PDF_submodule")
# Add the displacement/deformation functions
sys.path.append("./Displacement&Deformation_Functions")

# % Properties
cMapInd = 11; # Colourmap index (for later use with colour grading)
saveON = False; # Control flag for saving output
saveOutputDir = './Fig-EPS-PNG-output'  # The output directory

# % Initialize the mesh generator
m = MeshGenerator()
# % Setup the mesh
m.setup_mesh()
# % Plot the cube
tcube = m.translate(m.cube,np.array([-1, -1, -1]) * 0.5*0);
m.plot_elements(tcube);

# % Get the xyz mesh dimensions and properties
xyz = m.xyz;

# Define a function to update and plot the mesh with given displacement function and parameters
def run_wave_simulation(S, xyz_settings, colorF=[204/255, 204/255, 0], colorE=[100/255, 100/255, 100/255], face_alpha=1.0, plot_scale_factor=1.5, cMapInd=11, saveON=False, output_dir='./Fig-EPS-PNG-output', file_name_prefix='', mode='deformation'):
    """
    Run wave simulation with given settings and mesh configuration.
    
    Parameters:
    - S: Dictionary of settings for the displacement function.
    - xyz_settings: Mesh configuration settings.
    - colorF, colorE: Face and edge colors.
    - face_alpha: Transparency of faces.
    - plot_scale_factor: Scale factor for plotting.
    - cMapInd: Color map index for color grading.
    - saveON: Flag to save the output.
    - output_dir: Directory to save the output.
    - file_name_prefix: Prefix for the output file name.
    - mode: 'deformation' or 'translation' for applying transformations.
    """
    plt.close('all')  # Clear existing figures

    # Initialize the mesh generator and configure the mesh
    m = MeshGenerator()
    m.update_xyz_settings(xyz_settings)
    m.update_colors_and_plot_properties(colorF, colorE, face_alpha, plot_scale_factor)
    m.setup_mesh()

    # Get elements and apply displacements
    EE = m.get_elements()
    displacement_func = S['function']
    
    if mode == 'deformation':
        m.transform_elements(EE, displacement_func, S)
    elif mode == 'translation':
        for E in EE:
            displace_center_uxyz(E, S)

    # Apply color grading
    Es_original = m.get_elements()  # Assuming this captures the original state before transformation
    # Assuming E is an instance of CubeVoxel or similar,
    # and it has attributes `nodesR` and `id` you want to access.
    distM = np.max([np.linalg.norm(E.nodesR - Es_original[E.id].nodesR, axis=1) for E in EE])
    for E in EE:
        colour_grade(E, Es_original, distM, cMapInd)

    # Plot the mesh
    if cMapInd < 0:
        m.plot_mesh_fast(EE)
    else:
        m.plot_elements(EE)

    # Save the figure if saveON is True
    if saveON:
        file_name = f"{file_name_prefix}_{S['function'].__name__}_C_{cMapInd}.png"
        plt.savefig(Path(output_dir) / file_name)

    plt.draw()  # Display the plot
    plt.show()  # Display the plot
#end

## Example 1 : # Rayleigh Wave Simulation
S_rayleigh = {
    'function': rayleigh_wave_displacements,
    'c': 0.98,
    'cl': 0.85,
    'ct': 0.74,
    'A1': 3.4,
    'omega': 1,
    't': 0
}
xyz_settings_rayleigh = {
    'x': {'N': 50, 'gapP': 0.20, 'Lim': [-10, 10]},
    'y': {'N': 3,  'gapP': 0.20, 'Lim': [-1.5, 1.5]},
    'z': {'N': 20, 'gapP': 0.20, 'Lim': [-3, 0]}
}

run_wave_simulation(S_rayleigh, xyz_settings_rayleigh, file_name_prefix='RayleighWave', mode='deformation')
 


'''
## %% Example 1 = Make Mesh for RayleighWave

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
'''
