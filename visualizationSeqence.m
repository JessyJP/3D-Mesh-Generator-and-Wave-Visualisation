%% Sub functions 
function visualizationSeqence(m , S)
    global cMapInd saveON saveOutputDir exportCAD

    % Mesh properties
    % m.Pplot.axisEqual = true;
    % m.Pplot.View = 3;
    m.cube.colourF = [204 204 0]/255;
    m.cube.colourE = 70/100*[1 1 1];%[140 140 0]/255;
    m.cube.colourE = [100 100 100]/255;
    m.cube.faceAlpha = 100/100;

    % Plot Mesh
%     m.plotMeshFast(m.translate(m.elements,[0 0 0]*m.Pplot.scaleFactor));

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

    if exportCAD
        m.exportToFile(saveOutputDir+"/"+S.fileName,"STL")
        m.exportToFile(saveOutputDir+"/"+S.fileName,"OBJ")
    end
    disp("Done!")

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