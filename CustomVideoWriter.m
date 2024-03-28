classdef CustomVideoWriter < handle
    properties
        VideoObj % VideoWriter object
        SaveOutputDir % Directory where the video will be saved
        FPS % Frames per second
        Quality = 100; % Quality of the video
        f = 0;
    end
    
    methods
        function obj = CustomVideoWriter(videoName, saveOutputDir, fps)
            % Constructor for CustomVideoWriter
            obj.SaveOutputDir = saveOutputDir;
            obj.FPS = fps;
            videoFilename = fullfile(obj.SaveOutputDir, [videoName, '.mp4']);
            obj.VideoObj = VideoWriter(videoFilename, 'MPEG-4');
            obj.VideoObj.FrameRate = obj.FPS;
            obj.VideoObj.Quality = obj.Quality;
        end
        
        function open(obj)
            % Opens the video file for writing
            open(obj.VideoObj);
        end
        
        function addFrame(obj, mesh)
            % Adds a frame to the video
            % Ensure the frame is resized to match the desired output resolution if necessary
            drawnow;

            view(obj.f,obj.f)
        %     camproj('orthographic');
            camproj('perspective');
%             set(gca, 'CameraPosition', [0,0,0])
%             set(gca, 'CameraTarget', [0,0,1])
%             set(gca, 'CameraUpVector', [0,1,1])
        
            frame = getframe(mesh.fig); % Captures the figure as a frame
            resizedFrame = imresize(frame.cdata, [1080, 1920]); % Resize the frame
            writeVideo(obj.VideoObj, resizedFrame); % Writes the frame to the video file
    
%     close(fig); % Closes the figure to free up memory
            obj.f = obj.f + 1;
        end
        
        function close(obj)
            % Closes the video file
            close(obj.VideoObj);
        end
    end
end