# Temporal-segmentation
Long dissolve detection based in DBF

Preprocessing stage based in DBF

File name
main.m

Input
videoReader = VideoReader('input_video_name.mp4'); 
Output
videoWriter = VideoWriter('output_video_name.avi', 'Uncompressed AVI');

File name
temporalSegmentation2

Input
[GT_Vector, gtTransitions]= readGroundtruth('manual_anotations_GT_tool_based_format.txt',fps);
psnrs(1, :) = compareBhattHistsVideos([folder 'Preprocessed_DBF_video.avi'], [folder 'Non_preprocessed_video_Original.mp4'], gtTransitions);
Output
save('distanceResults');
