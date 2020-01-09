%make video from the frames (from input_folder); standard frame_rate = 30
function video_merge(input_folder, output) %output: name or full path; [AVI]
    working_dir = fullfile(pwd,input_folder);

    imageNames = dir(fullfile(working_dir,'*.jpg'));
    imageNames = {imageNames.name}';

    outputVideo = VideoWriter(output);
    outputVideo.FrameRate = 10; %!
    open(outputVideo);

    for ii = 1:length(imageNames)
       img = imread(fullfile(working_dir,imageNames{ii}));
       writeVideo(outputVideo,img);
    end

    close(outputVideo)
end




