%% calibration
C = imread('calibration.bmp');
figure; imshow(C);
imdistline
 
len = (32-10)/100;
px = 1777.6;
pxtom = len/px;

%% get particle centroids from images
fs = 10; % Hz
dt = 1/fs; % s
I_thres = 20; % particle intensity threshold
A_thres1 = 500; % min particle size threshold [px]
A_thres2 = 1500; % max particle size threshold [px]

d = dir; d(1:2) = [];
xp = zeros(length(d),1);
yp = xp;
figure; 
for i = 1:length(d)
    A = imread(d(i).name);
    imshow(A); hold on; 
    
    B = false(size(A));
    B(A < I_thres) = true;
    stats = regionprops('table',B,'centroid','area');
    idx = (stats.Area > A_thres1 & stats.Area < A_thres2);
    stats = stats(idx,:);
    if sum(idx) == 1
        xp(i) = stats.Centroid(:,1);
        yp(i) = stats.Centroid(:,2);
    else
        xp(i) = nan; yp(i) = nan;
    end
    
    plot(xp(i),yp(i),'r+','linewidth',2); pause(1/100); hold off
end

%% calculate rise velocity
r = sqrt(xp.^2 + (size(A,1) - yp).^2)*pxtom;
vel = gradient(r)/dt;

t = (0:length(d)-1)*dt;
figure; subplot(121); plot(t,r); subplot(122); plot(t,vel);

risevel = nanmean(vel(ceil(end/2):end-1));
fprintf('\nrise vel = %2.2f cm/s\n',risevel*100);
