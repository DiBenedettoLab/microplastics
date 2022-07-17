function [avar_k] = preprocess_MP(n)
% preprocess_MP

gdrive_path = 'C:\Users\ljbak\My Drive\';  %  'H:\My Drive\'; % 'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set
% n = 1;  % run number

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\run_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), ...
        run_params.ParticleType{n});

nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

load(sprintf('tracks_run%02d.mat',n))

% remove center bright spot
spot_x = [-.025 -.01];
spot_y = [-.3 -.28];

rm_idx = zeros(size(tracklength0));
for i = 1:length(tracklength0)
    idx = tracks0(:,5) == i;
    if sum(idx)
        if nanmean(tracks0(idx,1)) > spot_x(1) && nanmean(tracks0(idx,1)) < spot_x(2) && ...
                nanmean(tracks0(idx,2)) > spot_y(1) && nanmean(tracks0(idx,2)) < spot_y(2)
            rm_idx(i) = 1;
            tracks0(idx,:) = [];
            tracks0(tracks0(:,5)>i,5) = tracks0(tracks0(:,5)>i,5) - 1;
            i = i-1;
        end
    end
end
tracklength0(logical(rm_idx)) = [];

% repair broken tracks
searchrad = 10e-3;
[tracks0,tracklength0] = fix_tracks(tracks0,tracklength0,searchrad,1/run_params.imagingFreq_Hz(n),3);

% smooth tracks
sm_fn = sprintf('smtracks_run%02d.mat',n);
kernel = 5;% [3:2:21]; 
[smtracks, smtracklength, avar_k] = smooth_tracks(tracks0,kernel,1/run_params.imagingFreq_Hz(n));
if length(kernel) > 1
    figure; semilogy(kernel,avar_k,'k.'); 
    xlabel('kernel [frames]'), ylabel('var(a) [m^2/s^4]')
    k_opt = 5; k_opt_idx = logical(kernel >= k_opt);
    P = polyfit(kernel(k_opt_idx),log(avar_k(k_opt_idx)),1);
    hold on;
    plot(kernel,exp(kernel*P(1)+P(2)),'r-');
    plot(k_opt,avar_k(find(k_opt_idx,1,'first')),'ro','markersize',6);
    ylim([0 1.5])
end

save(sm_fn,'smtracks','smtracklength','kernel');

ntracks = length(smtracklength);

% get smoothed angles
if nonsphere
    kernel = 5;% 0;% 
    [smangles, smangles_cont] = get_smangles(tracks0,kernel,1/run_params.imagingFreq_Hz(n),run_params.ParticleType{n},run_params.Dp_m(n));
    save(sm_fn,'smangles','smangles_cont','-append');
end


%% preview tracks
% track lengths
figure; histogram(smtracklength,100)
xlabel('track length [frames]'); ylabel('count')

figure;
track_ids = round(linspace(1,ntracks,100));  % 1:30; %
c = jet(length(track_ids));
for i = 1:length(track_ids)
    idx = smtracks(:,5)==track_ids(i);
    c_idx = i; % round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));
    plot(smtracks(idx,1),smtracks(idx,2),'.','color',c(c_idx,:));
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')

