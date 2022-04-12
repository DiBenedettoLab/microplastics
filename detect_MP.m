% MP particle detection and coordinate mapping

% SETUP

clear
close all
warning off
addpath('H:\My Drive\MATLAB\fm-toolbox')
addpath('C:\Users\ljbak\My Drive\MATLAB\fm-toolbox')

% load experiment params
run_params = readtable('C:\Users\ljbak\My Drive\MP in OSBL\imaging expts\run_parameters_220315.ods');
cal_params = readtable('C:\Users\ljbak\My Drive\MP in OSBL\imaging expts\cal_parameters_220315.ods');
warning on

% run number
n = 1;
fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), ...
        run_params.ParticleType{n});

% make plots or not
plot_on = 1;

% save results or not
save_on = 0;

% image parameters 
cams = cell2mat(cal_params.Cam)';
dir_name = sprintf('run%i\\Cam%s\\', run_params.Run(n),cams(1));
imgset = dir([dir_name '*.tif']); 
img_nt = length(imgset);
img_shifts = floor([run_params.camA_offset_x(n), run_params.camA_offset_y(n); run_params.camC_offset_x(n), run_params.camC_offset_y(n); 
    run_params.camB_offset_x(n), run_params.camB_offset_y(n); run_params.camD_offset_x(n), run_params.camD_offset_y(n)]./2);

% adaptive binarization sensitivity
ABsens = run_params.adBinSensitivity(n); % 'Sensitivity' (range [0, 1]). High value thresholds more pixels as foregrnd, at risk of including some backgrnd pixels

% particle detection parameters
A_thres = [run_params.areaThreshold1_px(n), run_params.areaThreshold2_px(n)];

% set up figure
if plot_on
    Bfig = figure; 
    set(Bfig,'Units', 'Pixels','Position',[0.0010    0.0410    1.5360    0.7488]*1000); %[-1.9190   -0.1750    1.9200    0.9648]*1000);
    set(gcf,'color','w');
end



%% LOOP OVER CAMERAS

centers = cell(img_nt,length(cams)); % particle centroids [xp, yp]
angles = cell(img_nt,length(cams)); % particle orientation info [th_p, d_p]

for cam = 1:length(cams)   
    
    cam_left = cam <= 2; % 'true' for the left two cameras, 'false' for right two cameras
    
    % read image files from directory structure 
    dir_name = sprintf('run%i\\Cam%s\\', run_params.Run(n),cams(cam));
    imgset = dir([dir_name '*.tif']); 

    % load background and calibration images
    bkgd = cam_imread(sprintf('Cam%s-bkgd.tif',cams(cam)), cam_left);
    cal = cam_imread(sprintf('Cam%s-cal.tif',cams(cam)), cam_left);
    mask = cam_imread(sprintf('Cam%s-mask.tif',cams(cam)), cam_left);

    bkgd = double(bkgd);
    img_ix = size(bkgd,2); img_iy = size(bkgd,1);

    % crop to correct for camera shift
    bkgd_crop_rect = [abs(min([1,img_shifts(cam,1)])), max([1,img_shifts(cam,2)]), ...
        img_ix - max([0,img_shifts(cam,1)]), img_iy - abs(min([0,img_shifts(cam,2)]))];
    img_crop_rect = [max([1,img_shifts(cam,1)]), abs(min([1,img_shifts(cam,2)])), ...
        img_ix - abs(min([0,img_shifts(cam,1)])), img_iy - max([0,img_shifts(cam,2)])];
    
    bkgd = imcrop(bkgd, bkgd_crop_rect);
    cal = imcrop(cal, bkgd_crop_rect);
    mask = imcrop(mask, bkgd_crop_rect);

    img_ix = size(bkgd,2); img_iy = size(bkgd,1);



    % CALIBRATION: GET MAPPING FUNCTION FROM IMAGE TO WORLD COORDS
    
    cal = double(cal) - bkgd; % subtract background
    cal = uint8(cal - min(cal(:)));  % shift intensities so that all are positive 
    cal = 255 - cal;  % invert image
    
    %% binarize
    B = imbinarize(cal,'adaptive','Sensitivity',cal_params.calAdBinSens(cam)); % Adaptative binarization   
    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); % remove false positive pixels by eroding and dilating twice
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); 

    % number of rows and columns of dots in the calibration grid 
    n_rows_cal = cal_params.nRowsCal(cam); n_cols_cal = cal_params.nColsCal(cam);

    % mask top and edges (dots not part of the n_rows_cal x n_cols_cal grid)
    B = B.*logical(mask);
    if plot_on
        figure; subplot(161); pcolor_img(cal); title('bkgd sub and inverted')
        subplot(162); pcolor_img(B); title('masked and binarized')
    end

    % detect dots
    CC = bwconncomp(B);
    cal_dots = regionprops('table',CC,cal,'Centroid','Area','EquivDiameter');
    idx = cal_dots.Area > cal_params.calAreaThres1_px(cam) & cal_dots.Area < cal_params.calAreaThres2_px(cam);
%     idx = cal_dots.Area > 300 & cal_dots.Area < 5000;
    cal_dots = cal_dots(idx,:);    
    if plot_on
        subplot(163); pcolor_img(cal); hold on
        viscircles(cal_dots.Centroid,cal_dots.EquivDiameter/2); title('dots detected')
    end

    % known coords of the dots in physical space (world coords)
    W = generateCircleGridPoints([n_rows_cal, n_cols_cal], cal_params.spacing_m(cam), "PatternType","symmetric") + ...
        [cal_params.worldOffset_x(cam), cal_params.worldOffset_y(cam)]*cal_params.spacing_m(cam);

    % coords of the dots in image coords: sort dots into ascending order
    % from lower left corner of image by separating the dots by
    % x-coordinate into vertical bins   
    [~,top_row] = sort(cal_dots.Centroid(:,2),'descend'); 
    top_row = top_row(1:n_cols_cal);  % top row of dots
    bin_lim = cal_dots.Centroid(top_row,1); 
    bin_lim = sort(bin_lim,'ascend');
    bin_lim = [bin_lim - diff(bin_lim(1:2))/2; inf];  % bin limits
    
    I = zeros(n_rows_cal*n_cols_cal,2); 
    for j = 1:n_cols_cal    
        if plot_on; line([bin_lim(j) bin_lim(j)],[0 img_iy]); end  % plot bin limits
        cal_dots0 = cal_dots(cal_dots.Centroid(:,1) > bin_lim(j) & cal_dots.Centroid(:,1) < bin_lim(j+1),:);  % dots in a vertical bin
        [~,sort_idx] = sort(cal_dots0.Centroid(:,2),'ascend');  % sort by y-coord
        I(n_rows_cal*(j-1)+1 : n_rows_cal*j,:) = cal_dots0.Centroid(sort_idx,:);  % image points
    end
    if plot_on  % plot detected image points
        point_colors = jet(size(I,1));
        subplot(164); scatter(I(:,1),I(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('detected points')
        subplot(165); scatter(W(:,1),W(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('known coords')
    end

    % undistortion and rectification function
    rectify_quad = calibrate_camera(I,W,2);    
    
    % confirm that calibration dots are mapped to the correct world coords
    if plot_on
        imagePoints2 = rectify_quad(I); 
        subplot(166); scatter(imagePoints2(:,1),imagePoints2(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('corrected points')
    end
    


    %% LOOP OVER FRAMES
    i0 = 200; 
    nframes = 100; % img_nt; % 
    for i = i0:i0+nframes-1  
        if cam_left
            A = (imread([dir_name imgset(i).name]))'; % load image frame
        else
            A = rot90(imread([dir_name imgset(i).name]),2)'; 
        end
        A0 = imcrop(A,img_crop_rect);  % crop to correct for camera shift
        A0 = double(A0) - bkgd; % subtract background
        A0 = uint8(A0 - min(A0(:)));  % shift intensities so that all are positive 
        A0 = 255 - A0;  % invert image 
    
        % adaptive binarization
        B = imbinarize(A0,'adaptive','Sensitivity',ABsens); % Adaptative binarization   
        B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); % remove false positive pixels
        B = imerode(B,[0 1 0; 1 1 1; 0 1 0]);
        B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]);
        B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]);

        % free surface (develop this more)
        z_freesurf = inf;  % deepest extent of wave troughs [px]

        % FIND PARTICLES
        CC = bwconncomp(B);
        S1 = regionprops('table',CC,A0,'Centroid','Area','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','PixelIdxList');
        
        xp = []; yp = [];
        if ~isempty(S1)
            % remove based on area and proximity to free surface
            idx = S1.Area > A_thres(1) & S1.Area < A_thres(2) & S1.Centroid(:,2) < z_freesurf; % can also be a function of stats1.Centroid(:,1)
            S = S1(idx,:);
            CC.PixelIdxList = CC.PixelIdxList(idx);
    
            if ~isempty(S)
                % remove doubles (false shadows when the particle that shows up in addition to its shadow 
                % if particle is near the back wall) 
                dx_double = 300; % max horizontal distance between doubles
                dx_double_l = dx_double*~cam_left; % if particle is to the left of the shadow/right camera (px)
                dx_double_r = -dx_double*cam_left; % if particle is to the right of the shadow/left camera (px)
                dz_double = 20; % max vertical distance between doubles (px)
                Np = height(S);
                xp = S.Centroid(:,1);
                yp = S.Centroid(:,2);
                double_idx = ones(Np,1); 
                for j = 1:Np
                    distx = (xp(j) - xp);
                    disty = (yp(j) - yp);
                    for k = 1:Np
                        if j ~= k && distx(k) < dx_double_l && distx(k) > dx_double_r && disty(k) > -dz_double && disty(k) < dz_double 
                            double_idx(k) = 0; % flag doubles
                        end
                    end
                end            
                S = S(logical(double_idx),:);
                CC.PixelIdxList = CC.PixelIdxList(logical(double_idx));
    
                % particle centroids and orientations
                Np = height(S);
                xp = S.Centroid(:,1);
                yp = S.Centroid(:,2); 

                if strncmp(run_params.ParticleType{n},'d',1)
                    len1 = S.MajorAxisLength;
                    len2 = S.MinorAxisLength;
                    th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need (-S.Orientation) because image is flipped up-down)
                    d_p = len2; % apparent major axis length [px]
                elseif strncmp(run_params.ParticleType{n},'r',1)
                    len = sqrt(S.BoundingBox(:,3).^2 + S.BoundingBox(:,4).^2);
                    len = len - 0.25*S.MinorAxisLength;
                    th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need (-S.Orientation) because image is flipped up-down)
                    d_p = len; % apparent major axis length [px]
                end
            else
                Np = 0;
            end
        else
            Np = 0;
        end
    
        % accumulate centroids xp, yp in pixels and angles
        if Np > 0
            centers{i,cam} = [xp,yp];
            if strncmp(run_params.ParticleType{n},'r',1) || strncmp(run_params.ParticleType{n},'d',1)
                angles{i,cam} = [th_p,d_p];
            end
%             B2 = zeros(size(A)); B2(CC.PixelIdxList{1}) = 1;
        else
            centers{i,cam} = zeros(0,2);
            if strncmp(run_params.ParticleType{n},'r',1) || strncmp(run_params.ParticleType{n},'d',1)
                angles{i,cam} = zeros(0,2);
            end
        end
    
        % print progress
        if ~mod(i-1,50)
            fprintf([num2str(i-1) '/' num2str(img_nt) '\n'])
        end
    
        % MAKE PLOTS
        if plot_on
            figure(Bfig); clf;
            subplot(121); pcolor_img(B); % imadjust(A0,[200 255]/255,[150 255]/255)); hold on
            subplot(122); pcolor_img(A0); hold on
            if ~isempty(xp)
                % plot centroid
                plot(xp,yp,'r.','markersize',4'); 

                % plot major and minor axes
                if strncmp(run_params.ParticleType{n},'d',1)  
                    line([centers{i,cam}(:,1) - len1/2.*cos(th_p), centers{i,cam}(:,1) + len1/2.*cos(th_p)]', ...
                        [centers{i,cam}(:,2) - len1/2.*sin(th_p),centers{i,cam}(:,2) + len1/2.*sin(th_p)]', ...
                        'color','r','linewidth',0.75,'linestyle','-');
                    line([centers{i,cam}(:,1) - len2/2.*cos(th_p-pi/2), centers{i,cam}(:,1) + len2/2.*cos(th_p-pi/2)]', ...
                        [centers{i,cam}(:,2) - len2/2.*sin(th_p-pi/2),centers{i,cam}(:,2) + len2/2.*sin(th_p-pi/2)]', ...
                        'color','y','linewidth',0.75,'linestyle','-'); 
                elseif strncmp(run_params.ParticleType{n},'r',1)
                    line([centers{i,cam}(:,1) - len/2.*cos(th_p), centers{i,cam}(:,1) + len/2.*cos(th_p)]', ...
                        [centers{i,cam}(:,2) - len/2.*sin(th_p),centers{i,cam}(:,2) + len/2.*sin(th_p)]', ...
                        'color','r','linewidth',0.75,'linestyle','-'); 
                end
            end
%             set(gca,'XTick',[]); set(gca,'YTick',[]); 
            pause(1/100); 
            
%             % write to avi
%             if i == i0
%                 mname = 'test_disks_bkgd'; % movie file name 
%                 vid = VideoWriter(mname,'Grayscale AVI');
%                 vid.FrameRate = 10;
%                 open(vid)
%             end
%             writeVideo(vid,flipud(uint8(A0)));
%             if i == i0+nframes
%                 close(vid);
%             end
        
%             % write to gif
%             fn_gif = sprintf('C:\\Users\\ljbak\\My Drive\\MP in OSBL\\imaging expts\\%s-%i-%s.gif',run_params.ParticleType{n},run_params.WindSpeed_m_s(n),cams(cam));
%             fig_to_gif(fn_gif,0.1)
        end
   
    end

end


%% MERGE CAMERA VIEWS

% get array of coords of relevant points in px for each cam: endpoints of maj/min axes for disks,
% endpoints of minor axis for rods, centroid for nurdles/wax

% apply calibration: convert point coords in px into meters across all
% images

% convert points back into centroids and angles


%% SAVE CENTERS AND ANGLES
if save_on
    if strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1)
        save('centers.mat','centers','angles');
    else
        save('centers.mat','centers');
    end
end
