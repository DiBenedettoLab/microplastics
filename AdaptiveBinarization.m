%% -----------------BINARIZE IMAGES
% MANUAL INPUT=====================================================
fileLocat = 'C:\Users\ljbak\My Drive\MP in OSBL\basler\wasirf particle test 6 - nurdles'; % Location of the folder that contains the imgs
newFolder = 'Binarized'; % Name of new folder where binarized imgs will be saved
ABsens = 0.55; % 'Sensitivity' (range [0, 1]). High value thresholds more pixels as foregrnd, at risk including some backgrnd pixels
% -----------
% frameRate = 20; % frame rate for movie
% videoName = append('t4',newFolder,'Sen_',num2str(ABsens)); %name of video file
%==================================================================
mkdir(fileLocat, newFolder) % Make folder to store undistorted images
files = dir(fullfile(fileLocat, '*.tiff')); % Get all .tiff files
ImgNumber = length(files); % Give number of images

figure;
for id = 1:ImgNumber % Loop through each image
    I = imread(fullfile(fileLocat, files(id).name)); % Gets each image
    ABimage = imbinarize(I,'adaptive','Sensitivity',ABsens); % Adaptative binarization
    AB_uint8 = im2uint8(ABimage); % Converting back to uint8 so PIVlab can read the images
    imshow(AB_uint8);
    pause
%     imwrite(AB_uint8, fullfile(fileLocat, newFolder, files(id).name)); % Saves undistorted image in newFolder
end
MPEG4video(fullfile(fileLocat, newFolder), videoName, frameRate) % create the video writer
% clear I ABimage AB_uint8 fileLocat files ImgNumber newFolder id ABsens










return

%% JUST OTHER BINARIZATION PARAMETERS AND FUNCTIONS I TESTED OUT. 
%% test for different values of sensitivity in adaptive binarization
img = imread('/Users/juliochavez/Downloads/testImg/Undistorted/test_00000069.tiff');
ABimage1 = imbinarize(img,'adaptive','Sensitivity',0.05);
ABimage2 = imbinarize(img,'adaptive','Sensitivity',0.1);
ABimage3 = imbinarize(img,'adaptive','Sensitivity',0.2);
out = imtile({img, ABimage1, ABimage2, ABimage3}, 'Frames', 1:4, 'GridSize', [2 2]);
figure;
imshow(out);
% figure; imshowpair(imresize(img,0.5),imresize(ABimage,0.5),'montage');
% title('Original Image (left) vs. Corrected Image (right)');
clear img ABimage1 ABimage2 ABimage3 out
%% Adaptive histogram equalization
img2=adapthisteq(imgCropped); % apply adaptive hist equalization on previous output
img3=adapthisteq(img2); img4=adapthisteq(img3);
out = imtile({imgCropped, img2, img3, img4}, 'Frames', 1:4, 'GridSize', [2 2]);
figure;
imshow(out);
clear img imgCropped img2 img3 img4 out
% doesn't work because the image degenerates the more the adaptive
% histogram equalization algorithm is utilized.
% Pixel intensity comparison between foam and water reflexion doesn't work
% either because the intensity of the foam and the water reflexion are
% similar.
%% -----------------FUNCTIONS
function MPEG4video(fileLocation, videoName, frameRate)
files = dir(fullfile(fileLocation, '*.tiff')); % Get all .tiff files
ImgNumber = length(files); % Give number of images
video = VideoWriter(videoName, 'MPEG-4'); %create the video object
video.FrameRate = frameRate;
open(video); %open the file for writing
for ii=1:ImgNumber %where N is the number of images
  I = imread(fullfile(fileLocation, files(ii).name)); % Gets each image
  writeVideo(video,I); %write the image to file
end
close(video); %close the file
% clear video ii I files ImgNumber
end