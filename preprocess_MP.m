function preprocess_MP(n)
% preprocess_MP

gdrive_path = 'H:\My Drive\'; % C:\Users\ljbak\My Drive\';  %  'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220315';  % expt set
% n = 1;  % run number

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\run_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), ...
        run_params.ParticleType{n});

nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

load(sprintf('tracks_run%02d.mat',n))

% repair broken tracks
searchrad = 10e-3;
[tracks0,tracklength0] = fix_tracks(tracks0,tracklength0,searchrad,1/run_params.imagingFreq_Hz(n),3);

% smooth tracks
sm_fn = sprintf('smtracks_run%02d.mat',n);
kernel = 3;
[smtracks, smtracklength, avar_k] = smooth_tracks(tracks0,kernel,1/run_params.imagingFreq_Hz(n));
save(sm_fn,'smtracks','smtracklength','kernel');

ntracks = length(smtracklength);

% get smoothed angles
if nonsphere
    kernel = 3;
    [smangles, smangles_cont] = get_smangles(tracks0,kernel,1/run_params.imagingFreq_Hz(n),run_params.ParticleType{n},run_params.Dp_m(n));
    save(sm_fn,'smangles','smangles_cont','-append');
end


%% preview tracks
% track lengths
figure; histogram(smtracklength,100)
xlabel('track length [frames]'); ylabel('count')

figure;
track_ids = round(linspace(2,ntracks,100));  % 1:30; %
c = jet(length(track_ids));
for i = 1:length(track_ids)
    idx = smtracks(:,5)==track_ids(i);
    c_idx = i; % round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));
    plot(smtracks(idx,1),smtracks(idx,2),'.','color',c(c_idx,:));
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')

