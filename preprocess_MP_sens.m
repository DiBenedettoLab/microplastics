% preprocess_MP sensitivity

clear
close all

gdrive_path = 'C:\Users\ljbak\My Drive\';  %  'H:\My Drive\'; % 'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set
n = 4;  % run number

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\sens_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), ...
        run_params.ParticleType{n});

nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

load(sprintf('tracks_sens%02d.mat',n))

% % remove center bright spot
% spot_x = [-.025 -.01];
% spot_y = [-.3 -.28];
% 
% rm_idx = zeros(size(tracklength0));
% for i = 1:length(tracklength0)
%     idx = tracks0(:,5) == i;
%     if nanmean(tracks0(idx,1)) > spot_x(1) && nanmean(tracks0(idx,1)) < spot_x(2) && ...
%             nanmean(tracks0(idx,2)) > spot_y(1) && nanmean(tracks0(idx,2)) < spot_y(2)
%         rm_idx(i) = 1;
%         tracks0(idx,:) = [];
%     end
% end
% tracklength0(logical(rm_idx)) = [];

% % repair broken tracks
% searchrad = 10e-3;
% [tracks0,tracklength0] = fix_tracks(tracks0,tracklength0,searchrad,1/run_params.imagingFreq_Hz(n),3);

% % smooth tracks
% sm_fn = sprintf('smtracks_sens%02d.mat',n);
% kernel = 0;
% [smtracks, smtracklength, avar_k] = smooth_tracks(tracks0,kernel,1/run_params.imagingFreq_Hz(n));
% save(sm_fn,'smtracks','smtracklength','kernel');
% 
% ntracks = length(smtracklength);

% % get smoothed angles
% if nonsphere
%     kernel = 0;
%     [smangles, smangles_cont] = get_smangles(tracks0,kernel,1/run_params.imagingFreq_Hz(n),run_params.ParticleType{n},run_params.Dp_m(n));
%     save(sm_fn,'smangles','smangles_cont','-append');
% end

%% compute angles

Dp = run_params.Dp_m(n);
if strncmp(run_params.ParticleType{n},'r',1)
    % Rod angles (from raw tracks)
    th0r = tracks0(:,10);
    lr = tracks0(:,11);

    % subtract 1st-percentile l_r, scaled by rod angle
%     [N_lr,edges_lr] = histcounts(lr,200,'Normalization','cdf'); 
%     del_lr0 = mean(tracks0(:,12)); %5e-3; %edges_lr(find(N_lr>.99,1)+1) - Dp;  %min(lr); % [m]
%     del_lr = 1e-3;%mean(tracks0(:,12));
%     % del_lr = edges_lr(find(N_lr>.01,1)+1);  %min(lr); % [m]
% %     lr_shifted = max([(lr-del_lr.*(Dp-lr)/(Dp)), zeros(size(lr))],[],2);
%     lr_shifted = max([(lr - del_lr0), zeros(size(lr))],[],2);
%     lr_shifted = max([(lr_shifted - del_lr.*(Dp-lr_shifted)/(Dp)), zeros(size(lr_shifted))],[],2);
    del_lr0 = mean(tracks0(:,12));
    del_lr = 2.5e-3; % 1e-3;
    lr_shifted = lr - del_lr0;
    lr_shifted = lr_shifted - del_lr.*(Dp-lr_shifted)/Dp;
    
%     % plot lr histogram
%     [lr_pdf,lr_range] = pdf_var(lr*1000,100,0);
%     [lr_pdf2,lr_range2] = pdf_var(lr_shifted*1000,100,0);
%     figure; plot(lr_range,lr_pdf,'k.-',lr_range2,lr_pdf2,'r.-')
%     ylabel('PDF'); xlabel('$d$ [mm]'); ylim([0 1]); %legend('$d$','$d_{corr}$','location','nw'); xlim([0 3.5])
%     title('Rods')
%     goodplot([4 3.5])

    % compute angle cosines
    pxyz = [lr_shifted/(Dp).*cos(th0r), ... % p_x
    lr_shifted/(Dp).*sin(th0r), ... % p_y
    sqrt(1 - (lr_shifted/(Dp)).^2)]; % p_z
        
    rir = logical(~imag(pxyz(:,3)) & ~imag(pxyz(:,1))); % real idx, raw
    fprintf('%2.1f%% of raw p-hats real\n',sum(rir)/numel(rir)*100);   
    
elseif strncmp(run_params.ParticleType{n},'d',1)
    % Disk angles (from raw tracks)
    th0primer = tracks0(:,10);
    dr = tracks0(:,11);

    % subtract 1st-percentile d_r, scaled by disk angle
%     [N_dr,edges_dr] = histcounts(dr,200,'Normalization','cdf');
%     del_dr0 = mean(tracks0(:,12)) - Dp; 
%     del_dr = 1e-3; %edges_dr(find(N_dr>.01,1)+1);  %min(dr);
% %     del_dr = edges_dr(find(N_dr>.01,1)+1);  %min(dr);
%     dr_shifted = max([(dr - del_dr0), zeros(size(dr))],[],2);
%     dr_shifted = max([(dr_shifted - del_dr.*(Dp-dr_shifted)/(Dp)), zeros(size(dr_shifted))],[],2);
    del_dr0 = mean(tracks0(:,12)) - Dp; 
    del_dr = 0.5e-3; % 0;
    dr_shifted = dr - del_dr0;
    dr_shifted = dr_shifted - del_dr.*(Dp-dr_shifted)/Dp;
    
%     % plot dr histogram
%     [dr_pdf,dr_range] = pdf_var(dr*1000,100,0);
%     [dr_pdf2,dr_range2] = pdf_var(dr_shifted*1000,100,0);
%     figure; plot(dr_range,dr_pdf,'k.-',dr_range2,dr_pdf2,'r.-')
%     xlabel('$d$ [mm]'); legend('$d$','$d_{corr}$'); % xlim([0 2.5]);
%     ylim([0 1]); ylabel('PDF'); %set(gca,'YTickLabel',[]); %
%     title('Disks')
%     goodplot([4 3.5])
    
    % compute angle cosines
    pxyz = [sin(th0primer).*sqrt(1 - (dr_shifted/(Dp)).^2), ... % p_x
        cos(th0primer).*sqrt(1 - (dr_shifted/(Dp)).^2).*-sign(th0primer), ... % p_y
        dr_shifted/(Dp)]; % p_z
    
    rir = logical(~imag(pxyz(:,2)) & ~imag(pxyz(:,1))); % real idx, raw
    fprintf('%2.1f%% of raw p-hats real\n',sum(rir)/numel(rir)*100);
    
end

disp(pxyz)

% %% preview tracks
% % track lengths
% figure; histogram(smtracklength,100)
% xlabel('track length [frames]'); ylabel('count')
% 
% figure;
% track_ids = round(linspace(1,ntracks,100));  % 1:30; %
% c = jet(length(track_ids));
% for i = 1:length(track_ids)
%     idx = smtracks(:,5)==track_ids(i);
%     c_idx = i; % round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));
%     plot(smtracks(idx,1),smtracks(idx,2),'.','color',c(c_idx,:));
%     hold on
% end
% axis equal; axis([-.5 .5 -.45 .05]);
% xlabel('x [m]'); ylabel('y [m]')

