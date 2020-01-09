%output all frames into /output_folder
function video_explode(video, output_folder)
    working_dir = fullfile(pwd,output_folder); % ../output_folder
    mkdir(working_dir); %create the folder

    v = VideoReader(video);
    ii = 1;
    while hasFrame(v)
       img = readFrame(v);
       filename = [sprintf('%03d',ii) '.jpg'];
       fullname = fullfile(working_dir,filename);
       imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
       ii = ii+1;
    end
    disp("Video has been split into frames");
end




