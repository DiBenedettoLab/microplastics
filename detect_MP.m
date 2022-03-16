% MP PTV preprocessing
clear
close all
warning off
addpath('H:\My Drive\MATLAB\fm-toolbox')
addpath('C:\Users\ljbak\My Drive\MATLAB\fm-toolbox')
warning on

% read inputs from run_parameters.ods
plot_on = 1;
% load('tracks.mat')
% [~,tr_len_idx] = sort(tracklength0,'descend');
% tr_n = tr_len_idx(1);
% tr_idx = logical(tracks0(:,5) == tr_n | tracks0(:,5) == 5);
% plot_xy = [tracks0(tr_idx,1), tracks0(tr_idx,2)];
% xy_idx = tracks0(tr_idx,7);

% -- read image files from directory structure -- %
imgset = dir('*.tif'); 
img_nt = length(imgset);

%% background subtraction
  
% create background image (this should be done from a calibration image
% taken in still water)
bkgd = flipud(imread('../../CamA-bkgd2.tif')); 
img_ix = size(bkgd,2); img_iy = size(bkgd,1);
bkgd = int8(bkgd);


%% binarize (prelim processing)
centers = cell(img_nt,1); % particle centroids [xp, yp]
% angles = cell(img_nt,1); % particle orientation info [th_p, d_p]

% plot_xy = zeros(0,2);
% t_plot = zeros(0,1);

if plot_on
    Bfig = figure; 
    %set(Bfig,'Units', 'Pixels','Position',[2.6842   -0.0814    0.8080    0.7440]*1000);
    set(gcf,'color','w');
    pause
end

i0 = 100;
nframes = 150; % img_nt; % 
for i = i0:i0+nframes-1
    A = flipud(imread(imgset(i).name));
%     figure; imshow(A);
    A0 = int8(A) - bkgd;
    A0 = uint8(A0 - min(A0(:)));
    A0 = 255 - A0;
%     figure;pcolor_img(A0); keyboard

    %% adaptive binarization
    ABsens = 0.58; % 'Sensitivity' (range [0, 1]). High value thresholds more pixels as foregrnd, at risk including some backgrnd pixels
    B = imbinarize(A0,'adaptive','Sensitivity',ABsens); % Adaptative binarization

%     % non-adaptive binarization
%     thres = 235;
%     B = logical(A0 > thres);

    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); % remove tiny false positives
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]);
   
%     figure;pcolor_img(B)

    % detect particles
    % area
    A_thres1 = 300; A_thres2 = 2500; % nurdles
%     A_thres1 = 100; A_thres2 = 1800; % rods
%     A_thres1 = 200; A_thres2 = 2000; % disks
    CC = bwconncomp(B);
    stats1 = regionprops('table',CC,A0,'Centroid','Area','MeanIntensity','PixelIdxList');
    xp = []; yp = [];
    if ~isempty(stats1)
        % remove based on area and proximity to free surface
        idx = stats1.Area > A_thres1 & stats1.Area < A_thres2 & stats1.Centroid(:,2) < 1700;
        stats = stats1(idx,:);
        CC.PixelIdxList = CC.PixelIdxList(idx);

        if ~isempty(stats)
            % remove left and right false particles
            Np = height(stats);
            xp = stats.Centroid(:,1);
            yp = stats.Centroid(:,2);
            MI = stats.MeanIntensity;
            idx = ones(Np,1); 
            for j = 1:Np
                distx = (xp(j) - xp);
                disty = (yp(j) - yp);
                for k = 1:Np
                    if j ~= k && ... MI(k) < MI(j) && ...
                      distx(k) < 300 && distx(k) > 0 && disty(k) > -30 && disty(k) < 5 
                        idx(k) = 0;
                    end
                end
            end            
            stats = stats(logical(idx),:);
            CC.PixelIdxList = CC.PixelIdxList(logical(idx));

            Np = height(stats);
            xp = stats.Centroid(:,1);
            yp = stats.Centroid(:,2); 
        else
            Np = 0;
            xp = [];
            yp = [];
        end
    else
        Np = 0;
    end

    % found particles in frame
    if Np > 0

        % write centroids xp, yp in pixels
        centers{i} = [xp,yp];
%         angles{i} = [th_p,d_p];
        B2 = zeros(size(A)); B2(CC.PixelIdxList{1}) = 1;

%         if plot_on
% %             idx = logical(xp>1300 & yp>1000 & i<820);
% %             if any(idx)
% %             plot_xy = [plot_xy; xp(idx), yp(idx)];
%             plot_xy = [plot_xy; xp, yp];
%             t_plot = [t_plot;i];
% %             end
%         end
        
    % no particles in frame
    else
        centers{i} = zeros(0,2);
%         angles{i} = zeros(0,2);
    end

    % print progress
    if ~mod(i-1,200)
        fprintf([num2str(i-1) '/' num2str(img_nt) '\n'])
    end

    % plot
    if plot_on
        figure(Bfig); clf; pcolor_img(A0); hold on; % imadjust(A0,[200 255]/255,[150 255]/255)); hold on
        if ~isempty(xp)
            plot(xp,yp,'r.','markersize',4');
        end
%         if i >= xy_idx(1) && i <= xy_idx(end)
% %             plot(plot_xy(:,1),plot_xy(:,2),'r.','markersize',4)
%             plot(plot_xy(xy_idx<=(i-2),1),plot_xy(xy_idx<=(i-2),2),'r.','markersize',4)
%         end
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        pause(1/100); 
        
    %     % write to avi
    %     if i == i0
    %         mname = 'test_disks_bkgd'; % movie file name 
    %         vid = VideoWriter(mname,'Grayscale AVI');
    %         vid.FrameRate = 10;
    %         open(vid)
    %     end
    %     writeVideo(vid,flipud(uint8(A0)));
    %     if i == i0+nframes
    %         close(vid);
    %     end
    
%         % write to gif
%         cdata = print('-RGBImage','-r450');
%         frame.cdata = cdata; frame.colormap = [];
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if i == i0
%             imwrite(imind,cm,'disks_tracked.gif','gif', 'Loopcount',inf,'DelayTime',0.1);
%         else
%             imwrite(imind,cm,'disks_tracked.gif','gif','WriteMode','append','DelayTime',0.1);
%         end
    end
   
end

% % plotting arrays
% plot_xy_interp1 = interp1(t_plot,plot_xy(:,1),(i0:i)');
% plot_xy_interp2 = interp1(t_plot,plot_xy(:,2),(i0:i)');
% plot_xy = [plot_xy_interp1,plot_xy_interp2];


% % save centers
% save('centers.mat','centers'); %,'angles');
% 
% % track 
% kernel = 0;
% pxtom = 1;
% fs = 1;
% searchrad = 300;
% [tracks0,tracklength0] = ptvProcess2(centers,kernel,pxtom,fs,searchrad); %,angles);
% save('tracks.mat','tracks0','tracklength0','kernel')


%%
return























%% -- detect water surface/mask air layer -- %%


%% Get particle centroids 

centers = cell(img_nt,1); % particle centroids [xp, yp]
angles = cell(img_nt,1); % particle orientation info [th_p, d_p]

% intensity threshold for particle detection
if MP_shape == 'd'
    t1 = 150; t2 = 120; % detxn (5/29 2mm disks: t1 = 150; t2 = 120)
else
    t1 = 70; t2 = 60; % 70 60
end
part_thres = [linspace(t2,t1,img_ix/8),t1*ones(1,img_ix*6/8),linspace(t1,t2,img_ix/8)];   

fprintf('find particle centroids...\n')
thres = part_thres;
thres = repmat(thres,img_iy,1);
    
%  LOOP OVER FRAMES
for i = 1:img_nt
    A = imread(imgset(i).name);
    Ap = A - bkgd;

    % detect disks
    BW = zeros(size(Ap));
    BW(Ap > thres) = 1;

    % filter out small spurious objects
    CC = bwconncomp(logical(BW)); % connected bright regions
    if MP_shape == 'd'
        idx = find( cellfun(@numel,CC.PixelIdxList) < 0.6*2*Rp_px*thk/pxtom );
    else
        idx = find( cellfun(@numel,CC.PixelIdxList) < pi/4*(0.8*thk/pxtom)^2 ); 
    end
    for c = 1:length(idx)
        BW(CC.PixelIdxList{idx(c)}) = 0;
    end

    S = regionprops('table',logical(BW),Ap,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','PixelIdxList','MeanIntensity','Solidity');

    % reject partial particles, false positives 
    if ~isempty(S)
        if MP_shape == 'd'
            idx = ... %S.MinorAxisLength < 0.2*(2*Rp_px) & ...
                S.Centroid(:,1) < (img_ix - Rp_px) & S.Centroid(:,1) > Rp_px & S.Centroid(:,2) < (img_iy - Rp_px) & S.Centroid(:,2) > (y0+.2*Rp_px) & ...
                S.MajorAxisLength < 1.4*(2*Rp_px) & S.MajorAxisLength > 0.6*(2*Rp_px) & S.Solidity > .65 & ...
                S.MeanIntensity > 1.2*interp1(1:img_ix,thres(1,:),S.Centroid(:,1));
        else
            idx = S.Centroid(:,1) < (img_ix - Rp_px) & S.Centroid(:,1) > Rp_px & S.Centroid(:,2) < (img_iy - Rp_px) & ...
                    S.Area < 5*(2*Rp_px*thk/pxtom) & S.Area > 2*pi/4*(thk/pxtom)^2 & ...
                    S.MeanIntensity > 1.2*interp1(1:img_ix,thres(1,:),S.Centroid(:,1)) & S.Solidity > 0.75;
                len = sqrt(S.BoundingBox(:,3).^2 + S.BoundingBox(:,4).^2);
                idx = idx & (len < 2*1.2*Rp_px);
        end
        S2 = S(logical(idx),:);

        if MP_shape == 'd'
            len1 = S2.MajorAxisLength;
            len2 = S2.MinorAxisLength;
            xp = S2.Centroid(:,1); yp = S2.Centroid(:,2);
            th_p = -S2.Orientation(:)*pi/180; % in-plane orientation, radians (need (-S.Orientation) because image is flipped up-down)
            d_p = len2*pxtom; % apparent major axis length [m]
        else
            len = sqrt(S.BoundingBox(:,3).^2 + S.BoundingBox(:,4).^2);
            len = len - 0.25*S.MinorAxisLength;
            xp = S.Centroid(:,1); yp = S.Centroid(:,2);
            th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need (-S.Orientation) because image is flipped up-down)
            d_p = len*pxtom; % apparent major axis length [m]
        end

        Np = size(S,1);
    else
        Np = 0;
    end

    % found particles in frame
    if Np > 0

        % write centroids xp, yp in pixels
        centers{i} = [xp,yp];
        angles{i} = [th_p,d_p];
        Np = length(yp); 

    % no particles in frame
    else
        centers{i} = zeros(0,2);
        angles{i} = zeros(0,2);
    end

    % print progress
    if ~mod(i-1,1000)
        fprintf([num2str(i-1) '/' num2str(img_nt) '\n'])
    end

    if plot_on
        % plot detected particles
        figure(f)
        pcolor_img(A); title(num2str(i));
        if Np > 0    
            hold on; 
            if MP_shape == 'd'
                line([centers{i}(:,1) - len1/2.*cos(th_p), centers{i}(:,1) + len1/2.*cos(th_p)]', ...
                    [centers{i}(:,2) - len1/2.*sin(th_p),centers{i}(:,2) + len1/2.*sin(th_p)]', ...
                    'color','r','linewidth',1,'linestyle','-');
                line([centers{i}(:,1) - len2/2.*cos(th_p-pi/2), centers{i}(:,1) + len2/2.*cos(th_p-pi/2)]', ...
                    [centers{i}(:,2) - len2/2.*sin(th_p-pi/2),centers{i}(:,2) + len2/2.*sin(th_p-pi/2)]', ...
                    'color','b','linewidth',1,'linestyle','-'); 
    %            plot(centers{(i}(:,1), centers{i}(:,2),'r+','linewidth',1,'markersize',6);
            else
                line([centers{i}(:,1) - len/2.*cos(th_p), centers{i}(:,1) + len/2.*cos(th_p)]', ...
                        [centers{i}(:,2) - len/2.*sin(th_p),centers{i}(:,2) + len/2.*sin(th_p)]', ...
                        'color','r','linewidth',1,'linestyle','none','marker','+'); 
%                plot(centers{i}(:,1), centers{i}(:,2),'r+','linewidth',1,'markersize',6); 
            end
                hold off;

%             % write to gif
%             cdata = print('-RGBImage','-r600');
%             frame.cdata = cdata; frame.colormap = [];
%             im = frame2im(frame);
%             [imind,cm] = rgb2ind(im,256);
%             if i == i0
%                 imwrite(imind,cm,'nurdles_tracked.gif','gif', 'Loopcount',inf);
%             else
%                 imwrite(imind,cm,'nurdles_tracked.gif','gif','WriteMode','append');
%             end

%             % write to avi
%             if i == i0
%                 mname = 'test_nurdles'; % movie file name 
%                 vid = VideoWriter(mname,'Grayscale AVI');
%                 vid.FrameRate = 10;
%                 open(vid)
%             end
%             writeVideo(vid,A);            
    
            pause(1/1000) %keyboard %
        end
    end
end
    
% save centers
save('centers.mat','centers','angles');
