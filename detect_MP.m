% MP PTV preprocessing - disks
clear
close all
addpath('Z:\luci\MATLAB\pivio-matlab')
addpath('Z:\luci\MATLAB\PIV-PTV-toolbox')

load('input.mat')
plot_on = 0;

% -- read image files from directory structure instead of img -- %
img = imgRead([ifn '.img']);

%% background subtraction
  
% create background image
A = imgGetFrame(img,0);
N = 100;   % number of images to use
bkgd = zeros([size(A),2]);
bkgd(:,:,1) = A;
for i = round(linspace(1,img.it,N))
    bkgd(:,:,2) = imgGetFrame(img,i-1);
    bkgd(:,:,1) = min(bkgd,[],3);
end
bkgd = bkgd(:,:,1);

%%
if plot_on
    f = figure;
%     MP = get(0, 'MonitorPositions');
    set(f, 'Units', 'pixels');
    set(f, 'Position', [350 0 1534 1000]); %MP(1,:));
end

%% -- detect water surface/mask air layer -- %%


%% Get particle centroids 

centers = cell(img.it,1); % particle centroids [xp, yp]
angles = cell(img.it,1); % particle orientation info [th_p, d_p]

% intensity threshold for particle detection
if MP_shape == 'd'
    t1 = 150; t2 = 120; % detxn (5/29 2mm disks: t1 = 150; t2 = 120)
else
    t1 = 70; t2 = 60; % 70 60
end
part_thres = [linspace(t2,t1,img.ix/8),t1*ones(1,img.ix*6/8),linspace(t1,t2,img.ix/8)];   

fprintf('find particle centroids...\n')
thres = part_thres;
thres = repmat(thres,img.iy,1);
    
%  LOOP OVER FRAMES
for i = 1:img.it
    A = flipud(imgGetFrame(img,i-1));
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
                S.Centroid(:,1) < (img.ix - Rp_px) & S.Centroid(:,1) > Rp_px & S.Centroid(:,2) < (img.iy - Rp_px) & S.Centroid(:,2) > (y0+.2*Rp_px) & ...
                S.MajorAxisLength < 1.4*(2*Rp_px) & S.MajorAxisLength > 0.6*(2*Rp_px) & S.Solidity > .65 & ...
                S.MeanIntensity > 1.2*interp1(1:img.ix,thres(1,:),S.Centroid(:,1));
        else
            idx = S.Centroid(:,1) < (img.ix - Rp_px) & S.Centroid(:,1) > Rp_px & S.Centroid(:,2) < (img.iy - Rp_px) & ...
                    S.Area < 5*(2*Rp_px*thk/pxtom) & S.Area > 2*pi/4*(thk/pxtom)^2 & ...
                    S.MeanIntensity > 1.2*interp1(1:img.ix,thres(1,:),S.Centroid(:,1)) & S.Solidity > 0.75;
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
        fprintf([num2str(i-1) '/' num2str(img.it) '\n'])
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
%                 imwrite(imind,cm,'rods.gif','gif', 'Loopcount',inf);
%             else
%                 imwrite(imind,cm,'rods.gif','gif','WriteMode','append');
%             end
            pause(1/1000) %keyboard %
        end
    end
end
    
% save centers
save('centers.mat','centers','angles');
