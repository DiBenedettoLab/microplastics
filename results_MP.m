% process MP
% process microplastic particle positions/orientations

close all
clear

gdrive_path =  'C:\Users\ljbak\My Drive\';   %  'H:\My Drive\';  % 'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set
n = [15,13,1:6,16,14,7:12]; %  runs to include
runs_matched = [2 5 8 11 13 14]';  % runs with matching St (n3, r10, d7 runs)
n_triplet = find(any(n==runs_matched, 1));

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\run_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fs = run_params.imagingFreq_Hz(1);

nonsphere = zeros(size(n));
for i = 1:length(n)
    nonsphere(i) = strncmp(run_params.ParticleType{n(i)},'d',1) || strncmp(run_params.ParticleType{n(i)},'r',1);
end
nonsphere = logical(nonsphere);

smtracks = cell(size(n));
smtracklength = cell(size(n));
smangles = cell(size(n));
smangles_cont = cell(size(n));
for i = 1:(length(n))
    sm_tmp = load(sprintf('smtracks_run%02d.mat',n(i)));
    smtracks{i} = sm_tmp.smtracks;
    smtracklength{i} = sm_tmp.smtracklength;
    if nonsphere(i)
        smangles{i} = sm_tmp.smangles;
        smangles_cont{i} = sm_tmp.smangles_cont;
        imag_idx = logical(imag(smangles_cont{i}(:,1)));
        smangles{i}(imag_idx) = nan;
        smangles_cont{i}(imag_idx) = nan;
    end
end

% compute error bars (90% CI)
T_indep = 5;  % number of frames per independent realization
error_mean = @(q_var,q_N) 1.645*sqrt(q_var./(q_N/T_indep));
error_var = @(q_var,q_N) q_var.*(1-(q_N/T_indep-1)./chi2inv(0.05,q_N/T_indep-1));

% binning vars
Nbins = 30;%50;%12;   % number of bins in profiles
binedg = [-.4 -.05];
prof_zlim = [-.4 0];

Nbins_wide = 4;
z_bins_wide_edg = real(logspace(log10(binedg(1)),log10(binedg(2)),Nbins_wide+1)); % wide bin edges for autocorrelations and PDFs
z_bins_wide_cen = mean([z_bins_wide_edg(2:end); z_bins_wide_edg(1:end-1)],1);

Nt = zeros(size(n));
for i = 1:length(n); Nt(i) = max(smtracks{i}(:,7)); end

% plotting vars
ebar_gray = [.7 .7 .7]; 
ebar_red = [1 .6 .6]; 
ebar_blue = [.6 .6 1]; 
ebar_green = [0.6 0.9 0.6];
sym_green = [0 .5 0];

mk_col = {'k' 'k' 'k' 'r' 'r' 'b' 'b' 'b'};
mk_col2 = {ebar_gray ebar_gray ebar_gray ebar_red ebar_red ebar_blue ebar_blue ebar_blue};
mk_col = [mk_col mk_col2];
mk = {'.' '.' '.' '+' '+' 'o' 'o' 'o'};
mk = [mk mk];
mk_sz = [6 9 12 5 10 3 5 7];
mk_sz = [mk_sz mk_sz];
ls = {'-' '-' '-' '-' '-' '-' '-' '-'};
ls = [ls ls];
nones = {'none' 'none' 'none' 'none' 'none' 'none' 'none' 'none'};
nones = [nones nones];
eb_col = {ebar_gray ebar_gray ebar_gray ebar_red ebar_red ebar_blue ebar_blue ebar_blue};
eb_col = [eb_col eb_col];
lw = mk_sz./min(mk_sz)*1;
lstr = cellfun(@upper, run_params.ParticleType(n), 'UniformOutput', false)';  

% channel geometry & physical quantities
w_ROI = 0.91;
l_ROI = 0.98;
ROI_area = w_ROI*l_ROI; % horizontal area of ROI, channel width * ROI streamwise length
nu = 1e-6;
ut = 0.03;  % friction velocity
Hs = 0.05;  % significant wave height (eventually read from particle_flow_props spreadsheet?)
kdom = 1/.6;  % scaling (dominant wavenumber)
t_plus = 1;
AR = ones(size(n));
for i = 1:length(n)
    if strncmp(run_params.ParticleType{n(i)},'d',1)
        AR(i) = 0.8e-3/run_params.Dp_m(n(i));
    elseif strncmp(run_params.ParticleType{n(i)},'r',1)
        AR(i) = run_params.Dp_m(n(i))/1.75e-3;
    end
end

%% free surface elevation
if ~exist('phs-ztilde.mat','file')
    % wave phase
    
    % function to calculate phase at depth from phase at surface
    phs_depth = @(phs_surf, z) phs_surf;  % currently assuming phase at depth = phase at surface
    
    phs_tracks = cell(size(n));
    for i = 1:length(n)
        % load spatial phase data (from surface elevation)
        data = load(sprintf('C:/Users/ljbak/Documents/surf elev data/5 separated dewarped runs/run_%i_phase.mat',n(i)));
        x_phs = linspace(min(smtracks{i}(:,1)), max(smtracks{i}(:,1)), size(data.phs,1));
        
        % interpolate phase at each particle location
        phs_tracks{i} = zeros(length(smtracks{i}),1);
        for j = 1:max(smtracks{i}(:,7))  % loop over frames
            idx = logical(smtracks{i}(:,7) == j);   % get particles at frame j
    
            phs_tracks{i}(idx) = interp1(x_phs, data.phs(:,j), smtracks{i}(idx,1));  % interpolate phase at particle x coord
            phs_tracks{i}(idx) = phs_depth(phs_tracks{i}(idx), smtracks{i}(idx,2));  % correct phase for particle depth
        end
    end
    
    % free surface elev
    ztilde = cell(size(n));
    for i = 1:length(n)
        load(sprintf('C:/Users/ljbak/Documents/surf elev data/5 separated dewarped runs/run_%i_dewarped.mat',n(i)));
        ztilde{i} = zeros(length(smtracks{i}),1);
        for j = 1:size(smtracks{i},1)
            offset = interp1(surf_dewarp_x(:,smtracks{i}(j,7)), surf_dewarp_z(:,smtracks{i}(j,7)), smtracks{i}(j,1));
            ztilde{i}(j) = smtracks{i}(j,2) - offset;
        end
    end
    
    save('phs-ztilde.mat','phs_tracks','ztilde');
else
    load('phs-ztilde.mat')
end


%% particle spatial distribution
case_idx = 5+2;
N=hist3([smtracks{case_idx}(:,1),smtracks{case_idx}(:,2)],'Nbins',[100 50]);
figure;pcolor(N'); shading flat; axis equal tight


%% track length pdf
z_deep = -.1; % m
case_idx = 5+2;
deep_idx = smtracks{case_idx}(:,2) < z_deep;
uids_deep = unique(smtracks{case_idx}(deep_idx,5));

% uids_deep = [];
% for i = 1:length(smtracklength{case_idx})
%     idx = smtracks{case_idx}(:,5) == i;
%     if all(smtracks{case_idx}(idx,2) < z_deep)
%         uids_deep = [uids_deep; i];
%     end
% end

% % by frame count
% figure; histogram(smtracklength{case_idx}(uids),100)  % /run_params.imagingFreq_Hz(n(i)) [s]
% xlabel('Track length [frames]'); ylabel('Count')
% xlim([0 500]); goodplot([4 3])
% 
% length_long = 50;  % track length considered "long"
% fraction_long = sum(smtracklength{case_idx}(uids)>=length_long)/length(uids);  % fraction of tracks considered "long"
% fprintf('%2.1f%% of tracks deeper than %1.2f m are %2d frames or longer\n', fraction_long*100, -z_deep, length_long)

% by length in meters
% all tracks
trlen_m = zeros(size(smtracklength{case_idx}));
for i = 1:length(smtracklength{case_idx})
    idx = smtracks{case_idx}(:,5) == i; 
    tr_ds = sqrt(diff(smtracks{case_idx}(idx,1)).^2 + diff(smtracks{case_idx}(idx,2)).^2);  % track step lengths
    trlen_m(i) = sum(tr_ds); % track total length
end

% deep tracks
trlen_m = zeros(size(smtracklength{case_idx}));
for i = 1:length(smtracklength{case_idx})
    idx = smtracks{case_idx}(:,5) == i; 
    tr_ds = sqrt(diff(smtracks{case_idx}(idx,1)).^2 + diff(smtracks{case_idx}(idx,2)).^2);  % track step lengths
    trlen_m(i) = sum(tr_ds); % track total length
end

figure; histogram(trlen_m,60,'FaceColor',ebar_gray,'Normalization','pdf')  
xlabel('Track length [m]: all tracks'); ylabel('PDF'); xlim([0 1.5]); 
% annotation('textbox','string','(a)','position',[0.8 0.8 0.1 0.1],'edgecolor','none')
goodplot([4 4])

figure; histogram(trlen_m(uids_deep),60,'FaceColor',ebar_gray,'Normalization','pdf')  
xlabel(sprintf('Track length [m]: deeper than %1.1f m',z_deep)); ylabel('PDF'); xlim([0 1.5]); 
% annotation('textbox','string','(b)','position',[0.8 0.8 0.1 0.1],'edgecolor','none')
goodplot([4 4])

length_long = .2;  % track length considered "long"
fraction_long = sum(trlen_m(uids_deep)>=length_long)/length(uids_deep);  % fraction of tracks considered "long"
fprintf('%2.1f%% of tracks deeper than %1.2f m are %1.2f m or longer\n', fraction_long*100, -z_deep, length_long)



%% preview tracks by track number
figure;
track_ids = round(linspace(3,length(smtracklength{case_idx}),100));  % 1:100; %
c = turbo(length(track_ids)); % turbo(100); %
velmags = sqrt(smtracks{case_idx}(:,3).^2 + smtracks{case_idx}(:,4).^2);  %  smtracks{case_idx}(:,4); %  
velmags = velmags - min(velmags);
velmags = round(velmags*99/max(velmags)) + 1;
for i = 1:length(track_ids)
    idx = logical( smtracks{case_idx}(:,5)==track_ids(i) );
    c_idx = velmags(idx); % ceil(rand*length(track_ids));% round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));
    scatter(smtracks{case_idx}(idx,1),smtracks{case_idx}(idx,2),4,c(c_idx,:),'filled');
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')
hold off

goodplot([7 4.5])

%% preview tracks by frame number
% case_idx = 2+2;
% % load(sprintf('centers_run%02d.mat',n(case_idx)),'z_freesurf_inst_rect');
% load(sprintf('C:/Users/ljbak/Documents/surf elev data/5 separated dewarped runs/run_%i_dewarped.mat',n(case_idx)));
% %%
% figure; 
% set(gcf,'position',[369.8000  172.2000  817.6000  512.8000]); %[0.0010    0.0410    1.5360    0.7488]*1000)
% track_ids = round(linspace(2,length(smtracklength{1}),100));  % 1:30; %
% c = turbo(100); %turbo(length(track_ids));
% velmags = smtracks{case_idx}(:,4); %sqrt(smtracks{1}(:,3).^2 + smtracks{1}(:,4).^2); % 
% velmags = velmags - min(velmags);
% velmags = round(velmags*99/max(velmags)) + 1;
% stay_time = 60; % number of frames to keep on figure
% for i = 550%:800 %max(smtracks{1}(:,7))
%     idx1 = smtracks{case_idx}(:,7)==i; 
%     idx2 = smtracks{case_idx}(:,7)<i & smtracks{case_idx}(:,7)>(i-stay_time);
%     c_idx1 = velmags(idx1); %round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));  %   
%     c_idx2 = velmags(idx2); 
%     scatter(smtracks{case_idx}(idx1,1),smtracks{case_idx}(idx1,2),10,c(c_idx1,:),'filled'); hold on
%     scatter(smtracks{case_idx}(idx2,1),smtracks{case_idx}(idx2,2),2,c(c_idx2,:),'filled'); 
%     plot(-surf_dewarp_x(:,i),surf_dewarp_z(:,i),'k.','markersize',1); 
%     hold off
%     axis equal; axis([-.5 .5 -.45 .05]);
%     xlabel('x [m]'); ylabel('y [m]')
%     goodplot([7 5])
% 
%     pause(1/10)
% %     fig_to_gif('figs/tracks-r10-16.gif',0.1,'-r600')
% end


%% reconstructed tracks with orientation
case_idx = 6+2;
% uid = 1;

[~,track_idx_desc] = sort(smtracklength{case_idx},'descend');
uid = (track_idx_desc(40));  % particle ID % d5 50 r20 10 d7 202

figure;
plot_track_3d(smtracks{case_idx}, smangles{case_idx}, uid, run_params.ParticleType{n(case_idx)}, run_params.Dp_m(n(case_idx))/2);
set(gca,'ytick',''); ylabel('')
set(gca,'View',[28.9595 19.7060])
goodplot([6 4])


%% ------------------ EXIF ABOVE, RESULTS PAPER BELOW -------------------- %
% return































%% concentration

C = cell(size(n));
C_count = cell(size(n));
w_C = cell(size(n));
zprof = cell(size(n));
zprof_k = cell(size(n)); % z normalized by dom wavenumber
scaleflag = 0; % code not compatible with scaleflag=1!
for i = 1:length(n)
    % mean concentration
    [~, zprof{i}, C_count{i}] = condition_vars(ones(size(ztilde{i})),ztilde{i},Nbins,scaleflag,[-.3 0]);
%     [~, zprof{i}, C_count{i}] = condition_vars(ones(size(smtracks{i}(:,2))),smtracks{i}(:,2),Nbins,scaleflag,[-.4 0]); %binedg);
    C{i} = C_count{i}/(ROI_area*diff(zprof{i}(1:2))*Nt(i));
    zprof_k{i} = zprof{i}*kdom;
    
    % uncertainty
    C_tr = zeros(max(smtracks{i}(:,7)),Nbins);
    for j = 1:T_indep:Nt(i)
        idx = smtracks{i}(:,7) >= j & smtracks{i}(:,7) < j+T_indep;
        zp_j = smtracks{i}(idx,2);
        if ~isempty(zp_j)
            [~,~,C_tr(j,:)] = condition_vars(ones(size(zp_j)),zp_j,Nbins,scaleflag,binedg);
        end
    end
    C_tr = C_tr/(ROI_area*diff(zprof{i}(1:2))*T_indep);
    w_C{i} = error_mean(var(C_tr,1)',Nt(i));

    % scaling by uniform concentration
    if scaleflag
        normval = length(smtracks{i}(:,2))*real(diff(logspace(log10(binedg(1)),log10(binedg(2)),Nbins+1)'))/diff(binedg)/(ROI_area*diff(zprof{i}(1:2))*Nt(i));
    else
        normval = length(smtracks{i}(:,2))/Nbins/(ROI_area*diff(zprof{i}(1:2))*Nt(i));
    end
    Cnorm{i} = C{i}./normval;
    Cnorm{i}(Cnorm{i}==0) = 1e-4;
    w_C{i} = w_C{i}./normval;
end

%% plot concentration
figure;
l = compare_plots(Cnorm, zprof, mk_col, mk, mk_sz, ls,lw, ...
    w_C, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% set(gca,'yscale','log'); 
axis([0 4 prof_zlim])
xlabel('$C/C_0$'); ylabel('$z$ [m]'); legend(l,lstr,'location','se')
goodplot([6 5])

%% fit and theoretical profile
dz = 1*[0.07*ones(1,8), 0.055*ones(1,8)]; % 0.055; % zeros(1,8); %
Lm_fit = zeros(size(n));
C0 = zeros(size(n));
figure;
for i = length(n):-1:1
    idx = zprof{i} < -dz(i) & C_count{i}' > 0;  % remove wave layer and bottom of ROI
    P = polyfit((zprof{i}(idx)+dz(i)),log(Cnorm{i}(idx)),1);
    Lm_fit(i) = 1/P(1);
    C0(i) = exp(P(2));
    plot(log(Cnorm{i}(idx)),(zprof{i}(idx)+dz(i)),'-'); hold on; 
end

C0_count = zeros(size(n));
for i = 1:length(n)
    idx = C_count{i} > 0;
    P = polyfit((zprof{i}(idx)+dz(i)),log(C_count{i}(idx)),1);
    C0_count(i) = exp(P(2));
end

% % Lm vs Wb
% figure; plot(run_params.riseVel_m_s(n),Lm_fit,'+')
% xlabel('W_b [m/s]'); ylabel('L_m [m]')
% 
% % use fit Lm and measured rise vel to estimate diffusivity
% Lm_Wb = Lm_fit.*run_params.riseVel_m_s(n)';
% figure; plot(run_params.riseVel_m_s(n),Lm_Wb,'+')
% xlabel('W_b [m/s]'); ylabel('L_m*W_b [m^2/s]'); axis([0.015 0.06 2e-3 12e-3])

% theoretical Lm
A0 = 1.5*ut*0.4*Hs;
Lm_theor = A0./run_params.riseVel_m_s(n)';

% i = 1;
% z_exp = linspace(zprof{i}(1),zprof{i}(end),50);
% C_exp = C0(i)*exp(z_exp/Lm_fit(i));
% hold on; plot(C_exp,z_exp*kdom,'k--','linewidth',1)


% collapse conc profile
C_C0 = Cnorm;
w_C_C0 = w_C;
z_Lm = zprof;
for i = 1:length(n)
    C_C0{i} = Cnorm{i}/C0(i);
    w_C_C0{i} = w_C{i}/C0(i);
    z_Lm{i} = (zprof{i}+dz(i))/Lm_fit(i); % Lm_theor(i); %
end
figure; l = compare_plots(C_C0, z_Lm, mk_col, mk, mk_sz, ls, lw, ...
    w_C_C0, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
set(gca, 'Children', flipud(get(gca, 'Children')) )
hold on; l2 = plot(exp(-8:.1:0),-8:.1:0,'k--','linewidth',2);
ylim([-3 0]); xlim([0 2])
xlabel('$C/C_0$'); ylabel('$(z-\Delta z)/L_m$'); legend([l;l2],[lstr,{'Fit'}],'location','se')
goodplot([6 5])

fprintf('\nLm = %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f m\n',Lm_fit)
fprintf('C0 = %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f\n',C0)

% % check if developed 
% x_us = .3; x_ds = -x_us; us_idx = smtracks{1}(:,1) > x_us; ds_idx = smtracks{1}(:,1) < x_ds;
% [~, ~, C_us] = condition_vars(ones(size(smtracks{1}(us_idx,2))),smtracks{1}(us_idx,2),Nbins,scaleflag,binedg);
% [~, ~, C_ds] = condition_vars(ones(size(smtracks{1}(ds_idx,2))),smtracks{1}(ds_idx,2),Nbins,scaleflag,binedg);
% figure; plot(C_us, zprof{1}, 'b.', C_ds, zprof{1}, 'r.'); xlabel('C'); ylabel('z'); legend('upstream','downstream','location','se')


%% theoretical concentration profile, accounting for bottom wall

% theoretical C(z)
C_theor = cell(size(n));
H = 0.6;
z_theor = linspace(-H,0,30);
Lm_fit_test = 0.2;
for i = 1:length(n)
    C_C0_theor{i} = (1 + z_theor/(H-Lm_fit(i))).*exp(z_theor/Lm_fit(i));
%     if any(C_C0_theor{i}<0)  % if corrected profile goes negative,
% %     assumptions are incorrect, so revert to exponential decay (but: all
% %     are negative)
%         C_C0_theor{i} = exp(z_theor/Lm_fit(i)); 
%     end
end
% idx = 1:3;
% figure; compare_plots(C_C0_theor(idx), {z_theor z_theor z_theor}, mk_col(idx), {'none' 'none' 'none'} , mk_sz(idx), ls(idx),lw(idx));
% hold on; l = compare_plots(C_C0(idx), zprof(idx), mk_col(idx), mk(idx), mk_sz(idx), ls(idx), ...
%     w_C_C0(idx), num2cell(nan(size(n(idx)))), eb_col(idx), num2cell(nan(size(n(idx)))));
% legend(l,lstr,'location','se')
% hold on; l2 = plot(exp((-H:.01:0)/Lm_fit(idx(1))),-H:.01:0,'--','linewidth',2,'color','g');
% plot(exp((-H:.01:0)/Lm_fit(idx(2))),-H:.01:0,'--','linewidth',2,'color','g');
% if length(idx) > 2
%     plot(exp((-H:.01:0)/Lm_fit(idx(3))),-H:.01:0,'--','linewidth',2,'color','g');
% end
% % ylim([-H 0]); xlim([0 1.2])
% xlabel('$C/C_0$'); ylabel('$z$ [m]'); % legend([l;l2],[lstr,{'Fit'}],'location','se')
% goodplot([6 5])


%% vertical flux
% net flux
fluxz = cell(size(n));
for i = 1:length(n)
    flux_up = zeros(length(zprof{i}),1);
    flux_down = zeros(length(zprof{i}),1);
    for k = 1:length(zprof{i})
        [cross_up_idx, cross_down_idx] = flux_track(smtracks{i}, zprof{i}(k));
        flux_up(k) = sum(cross_up_idx);
        flux_down(k) = sum(cross_down_idx);
    end
    fluxz{i} = (flux_up - flux_down)/(ROI_area*Nt(i)/fs);
end

% concentration gradient and advective flux
Csmooth = cell(size(n));
dCdz = cell(size(n));
WC = cell(size(n));
for i = 1:length(n)
    Csmooth{i} = smooth(C{i},5); % smooth concentration profile
    dCdz{i} = gradient(Csmooth{i})./diff(zprof{i}(1:2));
    WC{i} = run_params.riseVel_m_s(n(i))*Csmooth{i};
end

% solve for diffusivity
eps_p_flux = cell(size(n));
for i = 1:length(n)
    eps_p_flux{i} = (WC{i} - fluxz{i})./dCdz{i};
end
A0_vec = 1.5*[1 2 3]*1e-2*0.4*5e-2;

% % plot flux balance components
% figure; l = compare_plots(WC,zprof_k,mk_col,mk,mk_sz,ls);
% xlabel('$W_b\langle C\rangle$'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([5 4])
% figure; l = compare_plots(dCdz,zprof_k,mk_col,mk,mk_sz,ls);
% xlabel('$d\langle C\rangle dz$'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([5 4])
% figure; l = compare_plots(fluxz,zprof_k,mk_col,mk,mk_sz,ls);
% xlabel('$\Phi$ [m$^{-2}$s$^{-1}$]'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([5 4])

% % plot diffusivity
% figure; l = compare_plots(eps_p_flux,zprof_k,mk_col,mk,mk_sz,ls);
% hold on; line([A0_vec(1) A0_vec(1)],binedg*kdom,'color','k','linestyle','--'); 
% line([A0_vec(2) A0_vec(2)],binedg*kdom,'color','k','linestyle','--'); 
% line([A0_vec(3) A0_vec(3)],binedg*kdom,'color','k','linestyle','--'); 
% xlim([-.01 .02]); xlabel('$A(z)$ [m$^{2}$s$^{-1}$]'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([4 4])



%% vertical excursions

delta_t_pdf = cell(length(n),1);
delta_t_rng = cell(length(n),1);
delta_zmax_pdf = cell(length(n),1);
delta_zmax_rng = cell(length(n),1);

Zc = -.06;  % crossing plane

for m = 1:length(n)
    delta_t = [];
    delta_zmax = [];

    [cross_up_idx, cross_down_idx] = flux_track(smtracks{m}, Zc);
    track_ids = unique(smtracks{m}(cross_up_idx | cross_down_idx,5));
    for j = 1:length(track_ids)
        t1 = find(cross_down_idx & smtracks{m}(:,5)==track_ids(j));
        t2 = find(cross_up_idx & smtracks{m}(:,5)==track_ids(j));
        if ~isempty(t1) && ~isempty(t2)
            t2 = t2(t2 > t1(1));
            for k = 1:length(t1)
                if length(t2) >= k
                    delta_t = [delta_t; t2(k) - t1(k)];
                    delta_zmax = [delta_zmax; min(smtracks{m}(t1(k):t2(k),2)) - Zc];
%                     keyboard
                end
            end
        end
    end
%         delta_zmax(delta_t>35) = [];
%         delta_t(delta_t>35) = [];
    delta_t = delta_t/fs;
    [delta_t_pdf{m}, delta_t_rng{m}] = pdf_var(delta_t,20,0); %1,[5/fs max(delta_t)]);
    [delta_zmax_pdf{m}, delta_zmax_rng{m}] = pdf_var(delta_zmax,20,0); %1,[min(delta_zmax) -1e-3]);
end

% % plot
% figure;
% compare_plots(delta_t_rng, delta_t_pdf, mk_col, mk, mk_sz, ls);
% xlabel('$\Delta t$ [s]'); ylabel('PDF'); set(gca,'yscale','log')
% goodplot([5 4])
% 
% figure;
% compare_plots(delta_zmax_rng, delta_zmax_pdf, mk_col, mk, mk_sz, ls);
% xlabel('$\Delta z_{max}$ [m]'); ylabel('PDF'); set(gca,'yscale','log')
% goodplot([5 4])



%% velocities

% pdfs
u_pdf = cell(size(n));
w_pdf = cell(size(n));
u_rng = cell(size(n));
w_rng = cell(size(n));

Npdf = 20; 

for i = 1:length(n)
    u_pdf{i} = zeros(Npdf,Nbins_wide);
    w_pdf{i} = zeros(Npdf,Nbins_wide);
    u_rng{i} = zeros(Npdf,Nbins_wide);
    w_rng{i} = zeros(Npdf,Nbins_wide);

    for j = 1:Nbins_wide
        idx = smtracks{i}(:,2) >= z_bins_wide_edg(j) & smtracks{i}(:,2) < z_bins_wide_edg(j+1);
        [u_pdf{i}(:,j), u_rng{i}(:,j)] = pdf_var(smtracks{i}(idx,3),Npdf,0);
        [w_pdf{i}(:,j), w_rng{i}(:,j)] = pdf_var(smtracks{i}(idx,4),Npdf,0);
    end
end

% plot
pdf_col = mat2cell(parula(Nbins_wide+1),ones(1,Nbins_wide+1),3);
pdf_col = pdf_col(1:end-1);
pdf_mk = {'.' '.' '.' '.' '.' '.'};
pdf_mk_sz = 10*ones(Nbins_wide,1);
pdf_ls = {'-' '-' '-' '-' '-' '-'};

% for i = 1:length(n)
%     figure;
%     subplot(121);
%     l = compare_plots(mat2cell(u_rng{i},Npdf,ones(1,Nbins_wide)), ...
%         mat2cell(u_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%     xlabel('$u_p$ [m/s]'); ylabel('PDF')
%     subplot(122);
%     compare_plots(mat2cell(w_rng{i},Npdf,ones(1,Nbins_wide)), ...
%         mat2cell(w_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%     xlabel('$w_p$ [m/s]')
%     sgtitle(run_params.ParticleType{i})
%     goodplot([5 4])
% end


% profiles
umean = cell(size(n));
wmean = cell(size(n));
uu = cell(size(n));
ww = cell(size(n));
uw = cell(size(n));
w_umean = cell(size(n));
w_wmean = cell(size(n));
w_uu = cell(size(n));
w_ww = cell(size(n));
w_uw = cell(size(n));
umean_total = zeros(size(n));
wmean_total = zeros(size(n));

zprof = cell(size(n));
zprof_k = cell(size(n)); % z normalized by dom wavenumber
scaleflag = 1;

for i = 1:length(n)
    [umean{i}, zprof{i}, N] = condition_vars(smtracks{i}(:,3),smtracks{i}(:,2),Nbins,scaleflag,binedg); 
    wmean{i} = condition_vars(smtracks{i}(:,4),smtracks{i}(:,2),Nbins,scaleflag,binedg); 
    zprof{i} = real(zprof{i});
    zprof_k{i} = zprof{i}*kdom;

    ufluct = smtracks{i}(:,3) - interp1(zprof{i},umean{i},smtracks{i}(:,2));
    wfluct = smtracks{i}(:,4) - interp1(zprof{i},wmean{i},smtracks{i}(:,2));
    uu{i} = condition_vars(ufluct.*ufluct,smtracks{i}(:,2),Nbins,scaleflag,binedg);
    ww{i} = condition_vars(wfluct.*wfluct,smtracks{i}(:,2),Nbins,scaleflag,binedg);
    uw{i} = condition_vars(ufluct.*wfluct,smtracks{i}(:,2),Nbins,scaleflag,binedg);

    w_umean{i} = error_mean(uu{i},N);
    w_wmean{i} = error_mean(ww{i},N);
    w_uu{i} = error_var(uu{i},N);
    w_ww{i} = error_var(ww{i},N);
    w_uw{i} = error_var(uw{i},N);

    umean_total(i) = mean(smtracks{i}(:,3),'omitnan');
    wmean_total(i) = mean(smtracks{i}(:,4),'omitnan');
end

%% plot velocity profiles
% figure;
% % subplot(121); 
% l = compare_plots(umean, zprof, mk_col, mk, mk_sz, ls, ...
%     w_umean, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% xlabel('$\langle u_p\rangle$ [m/s]'); ylabel('$z$ [m]')
% legend(l,lstr,'location','sw')
% 
% subplot(122); 
% compare_plots(wmean, zprof_k, mk_col, mk, mk_sz, ls, ...
%     w_wmean, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% xlabel('$\langle w_p\rangle$ [m/s]'); 
% goodplot([5 4])
% 
% figure;
% subplot(131); 
% l = compare_plots(uu, zprof_k, mk_col, mk, mk_sz, ls, ...
%     w_uu, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% xlabel('$\langle u_p''u_p''\rangle$ [m$^2$/s$^2$]'); ylabel('$zk_{dom}$')
% legend(l,lstr,'location','se')
% xlim([0 10]*1e-3)
% 
% subplot(132); 
% compare_plots(ww, zprof_k, mk_col, mk, mk_sz, ls, ...
%     w_ww, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% xlabel('$\langle w_p''w_p''\rangle$ [m$^2$/s$^2$]'); 
% xlim([0 10]*1e-3)
% 
% subplot(133); 
% compare_plots(uw, zprof_k, mk_col, mk, mk_sz, ls, ...
%     w_uw, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% xlabel('$\langle u_p''w_p''\rangle$ [m$^2$/s$^2$]'); 
% xlim([-5 5]*1e-3)
% goodplot([6 4])

% print total mean velocity
disp(umean_total)

figure; plot(run_params.riseVel_m_s(n), umean_total,'r.','markersize',10)
xlabel('$W_b$ [m/s]'); ylabel('$\langle \bar{u}\rangle$ [m/s]')
axis([0.01 0.06 0.03 0.07])
goodplot([5 4])


% % check if developed 
% umean_us = condition_vars(smtracks{1}(us_idx,3),smtracks{1}(us_idx,2),Nbins,scaleflag,binedg);
% wmean_us = condition_vars(smtracks{1}(us_idx,4),smtracks{1}(us_idx,2),Nbins,scaleflag,binedg);
% umean_ds = condition_vars(smtracks{1}(ds_idx,3),smtracks{1}(ds_idx,2),Nbins,scaleflag,binedg);
% wmean_ds = condition_vars(smtracks{1}(ds_idx,4),smtracks{1}(ds_idx,2),Nbins,scaleflag,binedg);
% figure; subplot(121); plot(umean_us, zprof{1}, 'b.', umean_ds, zprof{1}, 'r.'); xlabel('\langle u_p\rangle'); ylabel('z'); 
% subplot(122); plot(wmean_us, zprof{1}, 'b.', wmean_ds, zprof{1}, 'r.'); xlabel('\langle w_p\rangle'); ylabel('z'); legend('upstream','downstream','location','se')



%% orientations

smangles_nonan = smangles;


% orientation pdfs separated by depth
px_pdf = cell(size(n));
py_pdf = cell(size(n));
pz_pdf = cell(size(n));
pzs_pdf = cell(size(n));
px_rng = cell(size(n));
py_rng = cell(size(n));
pz_rng = cell(size(n));
pzs_rng = cell(size(n));

Npdf = 20; 
Npdfs = 1.5*Npdf+1;

for i = 1:length(n)
    if nonsphere(i)
        px_pdf{i} = zeros(Npdf,Nbins_wide);
        py_pdf{i} = zeros(Npdf,Nbins_wide);
        pz_pdf{i} = zeros(Npdf,Nbins_wide);
        pzs_pdf{i} = zeros(Npdfs,Nbins_wide);
        px_rng{i} = zeros(Npdf,Nbins_wide);
        py_rng{i} = zeros(Npdf,Nbins_wide);
        pz_rng{i} = zeros(Npdf,Nbins_wide);
        pzs_rng{i} = zeros(Npdfs,Nbins_wide);

        % fill in NaN's for rods
        if strncmp(run_params.ParticleType{n(i)},'r',1)
            nan_idx = isnan(smangles{i}(:,2));
            pz_nan = sin(smtracks{i}(:,10)).*nan_idx;
            smangles_nonan{i}(nan_idx,2) = pz_nan(nan_idx);
        end
    
        for j = 1:Nbins_wide 
            idx = smtracks{i}(:,2) >= z_bins_wide_edg(j) & smtracks{i}(:,2) < z_bins_wide_edg(j+1);
            [px_pdf{i}(:,j), px_rng{i}(:,j)] = pdf_var(abs(smangles{i}(idx,1)),Npdf,0,[0 1]);
            [py_pdf{i}(:,j), py_rng{i}(:,j)] = pdf_var(abs(smangles{i}(idx,3)),Npdf,0,[0 1]);
            [pz_pdf{i}(:,j), pz_rng{i}(:,j)] = pdf_var(abs(smangles_nonan{i}(idx,2)),Npdf,0,[0 1]);
            [pzs_pdf{i}(:,j), pzs_rng{i}(:,j)] = pdf_var(smangles_nonan{i}(idx,2),Npdfs,0,[-1 1]);
        end
    end
end

% theoretical PDF for random orientation distribution
pz_rng_rand = linspace(0,1,Npdf);
pz_pdf_rand = [];

% plot orientation pdfs
pdf_col = mat2cell(parula(Nbins_wide+1),ones(1,Nbins_wide+1),3);
pdf_col = pdf_col(1:end-1);
pdf_mk = {'.' '.' '.' '.' '.' '.'};
pdf_mk_sz = 10*ones(Nbins_wide,1);
pdf_ls = {'-' '-' '-' '-' '-' '-'};

for i = 1:length(n)
    if nonsphere(i)
        figure;
        subplot(141);
        l = compare_plots(mat2cell(px_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(px_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls,ones(size(pdf_mk_sz)));
        xlabel('$|p_x|$'); ylabel('PDF')
        subplot(142);
        compare_plots(mat2cell(py_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(py_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls,ones(size(pdf_mk_sz)));
        xlabel('$|p_y|$')
        subplot(143);
        compare_plots(mat2cell(pz_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(pz_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls,ones(size(pdf_mk_sz)));
        xlabel('$|p_z|$')
        subplot(144);
        compare_plots(mat2cell(pzs_rng{i},Npdfs,ones(1,Nbins_wide)), ...
            mat2cell(pzs_pdf{i},Npdfs,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls,ones(size(pdf_mk_sz)));
        xlabel('$p_z$')
        sgtitle(run_params.ParticleType{n(i)})
        goodplot([6 4])
    end
end

%% plot pz pdfs only
% pdf_col = mat2cell(parula(Nbins_wide),ones(1,Nbins_wide),3);
% pdf_mk = {'.' '.' '.' '.' '.' '.'};
% pdf_mk_sz = 10*ones(Nbins_wide,1);
% pdf_ls = {'-' '-' '-' '-' '-' '-'};

figure;
sp_count = 0;
for i = 1:length(n)
    if nonsphere(i)        
        sp_count = sp_count + 1;
        subplot(1,sum(nonsphere),sp_count);
        l = compare_plots(mat2cell(pz_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(pz_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
        xlabel('$|p_z|$')        
        title(upper(run_params.ParticleType{n(i)})) 
        ylim([0 5])
        if sp_count > 1
            set(gca,'yticklabel','')
        end
    end
end
goodplot([10 4])

%% one axis
% deep vs shallow
figure;
l = compare_plots(mat2cell([pz_rng{5}(:,[1,4]), pz_rng{7}(:,[1 4])],Npdf,ones(1,4)), ...
    mat2cell([pz_pdf{5}(:,[1 4]), pz_pdf{7}(:,[1 4])],Npdf,ones(1,4)), ...
    pdf_col([1 4 1 4]),{'+' '+' 'o' 'o'},[10 10 6 6],{'-' '-' '-' '-'},1.25*ones(1,5));
xlabel('$|p_z|$'); ylabel('PDF'); set(gca,'Yticklabel',[])
legend(l,{'Rods, deep', 'Rods, shallow', 'Disks, deep', 'Disks, shallow'},'location','nw')
ylim([0 7])
goodplot([5 5])

% all depths; compare to random orient distrib
figure;
l = compare_plots(mat2cell([pz_rng{4}(:,[1 end]), 1-pz_rng{7}(:,[1 end]), pz_rng{7}(:,1)],Npdf,ones(1,5)), ...
    mat2cell([pz_pdf{4}(:,[1 end]), pz_pdf{7}(:,[1 end]), ones(Npdf,1)],Npdf,ones(1,5)), ...
    [pdf_col([1 end],:); pdf_col([1 end],:); {[0 0 0]}],{'+' '+' 'o' 'o' 'none'},[10 10 6 6 1],{'-' '-' '-' '-' '-'},1.25*ones(1,5));
xlabel('$|p_z|$ or $1-|p_z|$'); ylabel('PDF'); %set(gca,'Yticklabel',[])
legend(l,{'Rods, deep', 'Rods, shallow', 'Disks, deep', 'Disks, shallow'},'location','ne')
ylim([0 5])
goodplot([5 5])


%% orient conditioned on asc and desc
if 0  % don't run this section
    asc_idx = cell(size(n));
    desc_idx = cell(size(n));
    for i = 1:length(n)
        asc_idx{i} = smtracks{i}(:,4) > 0;
        desc_idx{i} = smtracks{i}(:,4) < 0;
    end
    
    % orientation pdfs separated by depth  
    Npdf = 20; 
    Npdfs = 1.5*Npdf+1;
    
    for i = 1:length(n)
        if nonsphere(i)
            for j = 1:Nbins_wide
                idx = smtracks{i}(:,2) >= z_bins_wide_edg(j) & smtracks{i}(:,2) < z_bins_wide_edg(j+1);
                [px_pdf_asc{i}(:,j), px_rng_asc{i}(:,j)] = pdf_var(abs(smangles{i}(idx & asc_idx{i},1)),Npdf,0,[0 1]);
                [py_pdf_asc{i}(:,j), py_rng_asc{i}(:,j)] = pdf_var(abs(smangles{i}(idx & asc_idx{i},3)),Npdf,0,[0 1]);
                [pz_pdf_asc{i}(:,j), pz_rng_asc{i}(:,j)] = pdf_var(abs(smangles{i}(idx & asc_idx{i},2)),Npdf,0,[0 1]);
                [pzs_pdf_asc{i}(:,j), pzs_rng_asc{i}(:,j)] = pdf_var(smangles{i}(idx & asc_idx{i},2),Npdfs,0,[-1 1]);
    
                [px_pdf_desc{i}(:,j), px_rng_desc{i}(:,j)] = pdf_var(abs(smangles{i}(idx & desc_idx{i},1)),Npdf,0,[0 1]);
                [py_pdf_desc{i}(:,j), py_rng_desc{i}(:,j)] = pdf_var(abs(smangles{i}(idx & desc_idx{i},3)),Npdf,0,[0 1]);
                [pz_pdf_desc{i}(:,j), pz_rng_desc{i}(:,j)] = pdf_var(abs(smangles{i}(idx & desc_idx{i},2)),Npdf,0,[0 1]);
                [pzs_pdf_desc{i}(:,j), pzs_rng_desc{i}(:,j)] = pdf_var(smangles{i}(idx & desc_idx{i},2),Npdfs,0,[-1 1]);
            end
        end
    end
    
    for i = 1:length(n)
        if nonsphere(i)
            figure;
            subplot(141);
            l = compare_plots(mat2cell(px_rng_asc{i},Npdf,ones(1,Nbins_wide)), ...
                mat2cell(px_pdf_asc{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$|p_x|$'); ylabel('PDF')
            subplot(142);
            compare_plots(mat2cell(py_rng_asc{i},Npdf,ones(1,Nbins_wide)), ...
                mat2cell(py_pdf_asc{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$|p_y|$')
            subplot(143);
            compare_plots(mat2cell(pz_rng_asc{i},Npdf,ones(1,Nbins_wide)), ...
                mat2cell(pz_pdf_asc{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$|p_z|$')
            subplot(144);
            compare_plots(mat2cell(pzs_rng_asc{i},Npdfs,ones(1,Nbins_wide)), ...
                mat2cell(pzs_pdf_asc{i},Npdfs,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$p_z$')
            sgtitle([run_params.ParticleType{n(i)} ', ascending'])
            goodplot([6 4])
        end
    end
    
    for i = 1:length(n)
        if nonsphere(i)
            figure;
            subplot(141);
            l = compare_plots(mat2cell(px_rng_desc{i},Npdf,ones(1,Nbins_wide)), ...
                mat2cell(px_pdf_desc{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$|p_x|$'); ylabel('PDF')
            subplot(142);
            compare_plots(mat2cell(py_rng_desc{i},Npdf,ones(1,Nbins_wide)), ...
                mat2cell(py_pdf_desc{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$|p_y|$')
            subplot(143);
            compare_plots(mat2cell(pz_rng_desc{i},Npdf,ones(1,Nbins_wide)), ...
                mat2cell(pz_pdf_desc{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$|p_z|$')
            subplot(144);
            compare_plots(mat2cell(pzs_rng_desc{i},Npdfs,ones(1,Nbins_wide)), ...
                mat2cell(pzs_pdf_desc{i},Npdfs,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
            xlabel('$p_z$')
            sgtitle([run_params.ParticleType{n(i)} ', descending'])
            goodplot([6 4])
        end
    end
end


%% angular velocity (mean squared tumbling rate)
mstr = cell(size(n));
for i = 1:length(n)
    if nonsphere(i)
        sq_tumb_rate = sum(smangles{i}(:,4:6).^2,2);
        [mstr{i}, zprof2{i}] = condition_vars(sq_tumb_rate,smtracks{i}(:,2),20,0,prof_zlim);
    end
end

figure; 
compare_plots(mstr(nonsphere), zprof2(nonsphere), ...
    mk_col(nonsphere),nones(nonsphere),mk_sz(nonsphere),ls(nonsphere),lw(nonsphere));
ylabel('z [m]')
xlabel('$\langle\dot{p}\cdot\dot{p}\rangle$')
legend([lstr(nonsphere)],'location','northeastoutside')
goodplot([5 4])


%% depth & orientation = irradiation

% radiation as a function of depth 
lambda_list = linspace(.01,1.4,250); % [.05 .1 .2 .3 .4]; % .2; % light decay lengthscale [m]  %% can't compare with expt data if lambda < wave height due to lack of samples near surface
I0 = 1;  % surface light intensity
z_wave = 0;%-0.03; % wave layer depth

n_cases = length(n)*length(lambda_list);
Anormal = cell(n_cases,1);      % particle planar area normal to vertical, not normalized [m^2]
SA = zeros(n_cases,1);        % surface  area
Vp = zeros(n_cases,1);        % volume
Lm_lambda = zeros(n_cases,1);        % Lm_fit/lambda
irrad = cell(n_cases,1);      % irradiation [W]
irrad_depthonly = cell(n_cases,1); % irradiation on particle only considering depth
irrad_theor = zeros(n_cases,1);      % theoretical irradiation, depth-only [W?]
irrad_theor_rand = zeros(n_cases,1);      % theoretical irradiation, random orients [W?]
irrad_depthonly_extrap_mean = zeros(n_cases,1);
irrad_extrap_mean = zeros(n_cases,1);
w_irrad_depthonly_extrap_mean = zeros(n_cases,1);
w_irrad_extrap_mean = zeros(n_cases,1);
irrad_orientonly_rand_extrap_mean = zeros(n_cases,1);

smangles_nonan = smangles;


for j = 1:length(lambda_list)
    lambda = lambda_list(j);
    rad_level = @(z) I0*exp(z/lambda);   % radiation level as function of depth
    % rad_level = @(z) 1;

    depths = ztilde; % absolute or surface-relative particle depth
%     depths = cell(length(n),1);
%     for i = 1:length(n)
%         depths{i} = smtracks{i}(:,2);
%     end
    
    for i = 1:length(n)   
        idx = i + (j-1)*length(n);
        pz_rand = 0.5; % 3/8; % mean pz of random orientation distribution (use for disks)
        pxpy_rand = pi/4; % mean sqrt(1-pz^2) of random orientation distribution (use for rods)
        Lm_lambda(idx) = Lm_fit(i)/lambda;

        if nonsphere(i)  
            if strncmp(run_params.ParticleType{n(i)},'d',1)
                % disk
                thk = 0.8e-3;
                Ap = run_params.Dp_m(n(i))^2/4*pi; % particle nominal area
                A_eqsph = Ap*(3/2*AR(i))^(2/3);  % xsec area of equiv-volume sphere
                A_rand = Ap*pz_rand;  % mean normal area of randomly oriented particles
                SA(idx) = 2*pi*run_params.Dp_m(n(i))^2/4; % particle surface area
                Vp(idx) = pi*run_params.Dp_m(n(i))^2/4*thk; % particle volume
                SA_Vp(idx) = 2*(run_params.Dp_m(n(i)) + thk)/(run_params.Dp_m(n(i))*thk);  % SA to V ratio
                Anormal{idx} = abs(smangles{i}(:,2))*Ap; % normal area                
            else
                % rod
                thk = 1.75e-3;
                Ap = run_params.Dp_m(n(i))*thk; % particle nominal area
                A_eqsph = Ap*(9/16*1/AR(i))^(1/3);  % xsec area of equiv-volume sphere
                A_rand = Ap*pxpy_rand;  % mean normal area of randomly oriented rods
                SA(idx) = 2*pi*thk^2/4 + pi*thk*run_params.Dp_m(n(i)); % particle surface area
                Vp(idx) = pi*thk^2/4*run_params.Dp_m(n(i)); % particle volume
                SA_Vp(idx) = 2*(run_params.Dp_m(n(i)) + thk)/(run_params.Dp_m(n(i))*thk);  % SA to V ratio
                
                % replace NaN's with raw value of theta from images (which
                % is approx equal to pz for py~0, the case which results in NaN's)
%                 pxpy_interp = mean(sqrt(1 - (smangles{i}(smangles{i}(:,3)<.5,2).^2)),'omitnan'); 
                nan_idx = isnan(smangles{i}(:,2));
                pz_nan = sin(smtracks{i}(:,10)).*nan_idx;
                smangles_nonan{i}(nan_idx,2) = pz_nan(nan_idx);
                Anormal{idx} = sqrt(1 - (smangles_nonan{i}(:,2).^2))*Ap; % normal area 
            end
        else
            % sphere
            Ap = run_params.Dp_m(n(i))^2/4*pi;  % particle nominal area
            A_eqsph = Ap;   % xsec area of equiv-volume sphere
            A_rand = Ap;  % mean normal area of randomly oriented particles
            SA(idx) = 4*pi*run_params.Dp_m(n(i))^2/4; % particle surface area
            Vp(idx) = 4/3*pi*run_params.Dp_m(n(i))^3/8; % particle volume
            SA_Vp(idx) = 3/(run_params.Dp_m(n(i))/2);  % SA to V ratio
            Anormal{idx} = ones(length(depths{i}),1)*Ap;
        end
        
        % theoretical irradiance
        irrad_theor(idx) = I0*lambda/(lambda+Lm_fit(i))*exp(z_wave/lambda)/I0; % depth only (norm by surf irrad)
        irrad_theor_rand(idx) = irrad_theor(idx)*A_rand/(A_rand*I0); % random orientations (norm by rand particle at surface)
%         irrad_theor_depthlim(idx) = I0*lambda/(lambda+Lm_fit(i))*(1-exp(zprof{i}(1)/lambda)); % not correct, idk what's wrong
    
        % observed irradiance
        ROI_idx = smtracks{i}(:,2) > -.35; % % above bottom of ROI % true(size(depths{i}));  % smtracks{i}(:,2) < -dz(i);
        irrad_depthonly{idx} = rad_level(depths{i}(ROI_idx)); %
        irrad{idx} = Anormal{idx}(ROI_idx).*irrad_depthonly{idx};
        
        I_tot_extrap = C0_count(i)*Lm_fit(i)*irrad_theor(idx)*exp(zprof{i}(1)*(lambda+Lm_fit(i))/(lambda*Lm_fit(i)));
        N_tot_extrap = C0_count(i)*Lm_fit(i)*exp(zprof{i}(1)/Lm_fit(i));
        

        % depth-only irrad (norm by surface irrad)
        irrad_depthonly_extrap_mean(idx) = (sum(irrad_depthonly{idx},'omitnan') + I_tot_extrap) / ...
            (length(irrad_depthonly{idx}(~isnan(irrad_depthonly{idx}))) + N_tot_extrap) / I0;  

        % oriented irrad (normalized by rand orient particle at surface) 
        irrad_extrap_mean(idx) = (sum(irrad{idx},'omitnan') + I_tot_extrap*Ap) / ...
            (length(irrad{idx}(~isnan(irrad{idx}))) + N_tot_extrap) / (A_rand*I0);  
        
        % oriented irrad (not normalized) 
        irrad_extrap_mean_abs(idx) = (sum(irrad{idx},'omitnan') + I_tot_extrap*Ap) / ...
            (length(irrad{idx}(~isnan(irrad{idx}))) + N_tot_extrap) / I0;  

        % depth-only irrad assuming random-orient 
        irrad_depthrand_extrap_mean(idx) = (sum(irrad_depthonly{idx},'omitnan') + I_tot_extrap)*A_rand / ... 
            (length(irrad_depthonly{idx}(~isnan(irrad_depthonly{idx}))) + N_tot_extrap) / I0;  
        
        % depth-only irrad assuming perfect pref orient 
        irrad_depthpref_extrap_mean(idx) = (sum(irrad_depthonly{idx},'omitnan') + I_tot_extrap)*Ap / ... 
            (length(irrad_depthonly{idx}(~isnan(irrad_depthonly{idx}))) + N_tot_extrap) / I0; 

        % orient-only irrad (oriented irrad normalized by depth-only random irrad)
        irrad_orientonly_rand_extrap_mean(idx) = irrad_extrap_mean_abs(idx) ./ irrad_depthrand_extrap_mean(idx);

        % orient-only irrad (oriented irrad normalized by depth-only perf pref irrad)
        irrad_orientonly_pref_extrap_mean(idx) = irrad_extrap_mean_abs(idx) ./ irrad_depthpref_extrap_mean(idx);
        
%         ( (sum(irrad{idx},'omitnan') + I_tot_extrap*Ap) / ...
%             (length(irrad{idx}(~isnan(irrad_depthonly{idx}))) + N_tot_extrap) ) ./ ...
%             ( (sum(irrad_depthonly{idx},'omitnan') + I_tot_extrap)*A_rand / ...
%             (length(irrad_depthonly{idx}(~isnan(irrad_depthonly{idx}))) + N_tot_extrap) );

        % error estimates [TODO]
        w_irrad_depthonly_extrap_mean(idx) = nan;
        w_irrad_extrap_mean(idx) = nan;
    end
end


% % irradiation vs Lm considering depth only, no orientation
% figure;
% compare_plots(num2cell(Lm_lambda), num2cell(irrad_depthonly_extrap_mean), ...
%     repmat(mk_col,1,length(lambda_list)),repmat(mk,1,length(lambda_list)),repmat(mk_sz,1,length(lambda_list)),repmat(nones,1,length(lambda_list))); hold on; 
% xlabel('$L_m/\lambda$ [m]'); ylabel('$I/I_0$')
% legend(lstr,'location','northeastoutside')
% % P = polyfit(log(Lm_fit/lambda),log(irrad_depthonly_extrap_mean),1);
% % Lm_vec = 0.4:0.01:2;
% % plot(Lm_vec, exp(P(2))*Lm_vec.^P(1));
% plot(Lm_lambda, irrad_theor, '*','color',ebar_gray);
% title('Extrapolated: depth only')
% goodplot([6 4])
% 
% % irradiation vs Lm accounting for orientation, normalized by random pz distribution
% figure;
% compare_plots(num2cell(Lm_lambda), num2cell(irrad_extrap_mean), ...
%     repmat(mk_col,1,length(lambda_list)),repmat(mk,1,length(lambda_list)),repmat(mk_sz,1,length(lambda_list)),repmat(nones,1,length(lambda_list))); hold on; 
% xlabel('$L_m/\lambda$ [m]'); ylabel('$I/I_0$')
% legend(lstr,'location','northeastoutside')
% % P = polyfit(log(Lm_fit/lambda),log(irrad_depthonly_extrap_mean),1);
% % Lm_vec = 0.4:0.01:2;
% % plot(Lm_vec, exp(P(2))*Lm_vec.^P(1));
% plot(Lm_lambda, irrad_theor, '*','color',ebar_gray);
% title('Extrapolated: orient + depth')
% goodplot([6 4])

% line plots
% irradiation vs Lm considering depth only, no orientation
figure;
compare_plots(num2cell(reshape(Lm_lambda,length(n),[]),2), num2cell(reshape(irrad_depthonly_extrap_mean,length(n),[]),2), ...
    mk_col,nones,mk_sz,ls,lw); hold on; 
xlabel('$L_m/\lambda$ [m]'); ylabel('$\langle I\rangle/I_0$')
legend(lstr,'location','northeastoutside')
% P = polyfit(log(Lm_fit/lambda),log(irrad_depthonly_extrap_mean),1);
% Lm_vec = 0.4:0.01:2;
% plot(Lm_vec, exp(P(2))*Lm_vec.^P(1));
plot(Lm_lambda, irrad_theor, '.','color',ebar_gray);
title('Extrapolated: depth only')
axis([0 6 0 1.5])
goodplot([6 4])

% irradiation vs Lm accounting for orientation, normalized by random pz distribution
figure;
compare_plots(num2cell(reshape(Lm_lambda,length(n),[]),2), num2cell(reshape(irrad_extrap_mean,length(n),[]),2), ...
    mk_col,nones,mk_sz,ls,lw); hold on; 
xlabel('$L_m/\lambda$ [m]'); ylabel('$\langle I\rangle/I_0$')
legend(lstr,'location','northeastoutside')
% P = polyfit(log(Lm_fit/lambda),log(irrad_depthonly_extrap_mean),1);
% Lm_vec = 0.4:0.01:2;
% plot(Lm_vec, exp(P(2))*Lm_vec.^P(1));
plot(Lm_lambda, irrad_theor_rand, '.','color',ebar_gray);
title('Extrapolated: orient + depth')
axis([0 6 0 1.5])
goodplot([6 4])

% irradiation vs Lm accounting for orientation, normalized by depth-only rand orient irrad
Lm_lambda_cell = num2cell( reshape(Lm_lambda,length(n),[]), 2 );
irrad_cell = num2cell(reshape(irrad_orientonly_rand_extrap_mean,length(n),[]),2);
figure;
compare_plots(Lm_lambda_cell(nonsphere), irrad_cell(nonsphere), ...
    mk_col(nonsphere),nones(nonsphere),mk_sz(nonsphere),ls(nonsphere),lw(nonsphere)); hold on; 
xlabel('$L_m/\lambda$'); ylabel('$\langle I\rangle/I_{rand}$')
line([0 6],[1 1],'color','k');
legend([lstr(nonsphere) {'Random'}],'location','northeastoutside')
% P = polyfit(log(Lm_fit/lambda),log(irrad_depthonly_extrap_mean),1);
% Lm_vec = 0.4:0.01:2;
% plot(Lm_vec, exp(P(2))*Lm_vec.^P(1));
% plot(Lm_lambda, irrad_theor_rand, '.','color',ebar_gray);
% title('Extrapolated: orient only (rand)')
axis([0 6 0.9 1.4])
goodplot([6 4])

% irradiation vs Lm accounting for orientation, normalized by depth-only perf pref orient irrad
irrad_cell = num2cell(reshape(irrad_orientonly_pref_extrap_mean,length(n),[]),2);
figure;
compare_plots(Lm_lambda_cell(nonsphere), irrad_cell(nonsphere), ...
    mk_col(nonsphere),nones(nonsphere),mk_sz(nonsphere),ls(nonsphere),lw(nonsphere)); hold on; 
xlabel('$L_m/\lambda$ [m]'); ylabel('$\langle I\rangle/I_0$')
line([0 6],[1 1],'color','k');
legend([lstr(nonsphere) {'Random'}],'location','northeastoutside')
% P = polyfit(log(Lm_fit/lambda),log(irrad_depthonly_extrap_mean),1);
% Lm_vec = 0.4:0.01:2;
% plot(Lm_vec, exp(P(2))*Lm_vec.^P(1));
% plot(Lm_lambda, irrad_theor_rand, '.','color',ebar_gray);
title('Extrapolated: orient only (pref)')
axis([0 6 0 1.5])
goodplot([6 4])




%% random walk model
return;
addpath([gdrive_path 'MATLAB\model'])

irrad_rw = zeros(size(n));
nu_ts = Lm_fit.*run_params.riseVel_m_s(n)';
for i = 1:length(n)
    irrad_rw(i) = random_walk_irrad(run_params.riseVel_m_s(i), nu_ts(i), lambda, Lm_fit(i), H);
end

figure;
compare_plots(num2cell(Lm_fit/lambda), num2cell(irrad_depthonly_extrap_mean),mk_col,mk,mk_sz,nones); hold on; 
xlabel('$L_m/\lambda$ [m]'); ylabel('$I/I_0$')
goodplot([5 4])
plot(Lm_fit/lambda, irrad_theor, '*','color',ebar_gray);

compare_plots(num2cell(Lm_fit/lambda), num2cell(irrad_rw),mk_col,mk,mk_sz,nones); hold on; 
xlabel('$L_m/\lambda$ [m]'); ylabel('$I/I_0$')
goodplot([5 4])
plot(Lm_fit/lambda, irrad_theor, '*','color',ebar_gray);
title('Random walk')



%% phase averaging

% % condition orientation on wave phase 
% % PDFs
% Nbins_phs = 7;
% phs_bins_edg = linspace(0,2*pi,Nbins_phs+1);
% 
% px_pdf_phs = cell(size(n));
% py_pdf_phs = cell(size(n));
% pz_pdf_phs = cell(size(n));
% pzs_pdf_phs = cell(size(n));
% px_rng_phs = cell(size(n));
% py_rng_phs = cell(size(n));
% pz_rng_phs = cell(size(n));
% pzs_rng_phs = cell(size(n));
% 
% Npdf = 20; 
% Npdfs = 1.5*Npdf+1;
% 
% for i = 1:length(n)
%     if nonsphere(i)
%         px_pdf_phs{i} = zeros(Npdf,Nbins_phs);
%         py_pdf_phs{i} = zeros(Npdf,Nbins_phs);
%         pz_pdf_phs{i} = zeros(Npdf,Nbins_phs);
%         pzs_pdf_phs{i} = zeros(Npdfs,Nbins_phs);
%         px_rng_phs{i} = zeros(Npdf,Nbins_phs);
%         py_rng_phs{i} = zeros(Npdf,Nbins_phs);
%         pz_rng_phs{i} = zeros(Npdf,Nbins_phs);
%         pzs_rng_phs{i} = zeros(Npdfs,Nbins_phs);
%     
%         for j = 1:Nbins_phs
%             idx = phs_tracks{i} >= phs_bins_edg(j) & phs_tracks{i} < phs_bins_edg(j+1);% & smtracks{i}(:,2) < -.2;
%             [px_pdf_phs{i}(:,j), px_rng_phs{i}(:,j)] = pdf_var(abs(smangles{i}(idx,1)),Npdf,0,[0 1]);
%             [py_pdf_phs{i}(:,j), py_rng_phs{i}(:,j)] = pdf_var(abs(smangles{i}(idx,3)),Npdf,0,[0 1]);
%             [pz_pdf_phs{i}(:,j), pz_rng_phs{i}(:,j)] = pdf_var(abs(smangles{i}(idx,2)),Npdf,0,[0 1]);
%             [pzs_pdf_phs{i}(:,j), pzs_rng_phs{i}(:,j)] = pdf_var(smangles{i}(idx,2),Npdfs,0,[-1 1]);
%         end
%     end
% end
% 
% % plot params
% pdf_col = mat2cell(parula(Nbins_phs+1),ones(1,Nbins_phs+1),3);
% pdf_col = pdf_col(1:end-1);
% pdf_mk = {'.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.'};
% pdf_mk_sz = 10*ones(Nbins_phs,1);
% pdf_ls = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-'};
% 
% % % plot orientation pdfs 
% % for i = 1:length(n)
% %     if nonsphere(i)
% %         figure;
% %         subplot(141);
% %         l = compare_plots(mat2cell(px_rng_phs{i},Npdf,ones(1,Nbins_phs)), ...
% %             mat2cell(px_pdf_phs{i},Npdf,ones(1,Nbins_phs)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
% %         xlabel('$|p_x|$'); ylabel('PDF')
% %         subplot(142);
% %         compare_plots(mat2cell(py_rng_phs{i},Npdf,ones(1,Nbins_phs)), ...
% %             mat2cell(py_pdf_phs{i},Npdf,ones(1,Nbins_phs)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
% %         xlabel('$|p_y|$')
% %         subplot(143);
% %         compare_plots(mat2cell(pz_rng_phs{i},Npdf,ones(1,Nbins_phs)), ...
% %             mat2cell(pz_pdf_phs{i},Npdf,ones(1,Nbins_phs)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
% %         xlabel('$|p_z|$')
% %         subplot(144);
% %         compare_plots(mat2cell(pzs_rng_phs{i},Npdfs,ones(1,Nbins_phs)), ...
% %             mat2cell(pzs_pdf_phs{i},Npdfs,ones(1,Nbins_phs)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
% %         xlabel('$p_z$')
% %         sgtitle(run_params.ParticleType{n(i)})
% %         goodplot([6 4])
% %     end
% % end
% 
% % plot pz pdfs only
% figure;
% sp_count = 0;
% for i = 1:length(n)
%     if nonsphere(i)        
%         sp_count = sp_count + 1;
%         subplot(1,sum(nonsphere),sp_count);
%         l = compare_plots(mat2cell(pz_rng_phs{i},Npdf,ones(1,Nbins_phs)), ...
%             mat2cell(pz_pdf_phs{i},Npdf,ones(1,Nbins_phs)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%         xlabel('$|p_z|$')        
%         title(upper(run_params.ParticleType{n(i)})) 
%         ylim([0 5])
%         if sp_count > 1
%             set(gca,'yticklabel','')
%         else
%             ylabel('PDF')
%         end
%     end
% end
% goodplot([10 4])


%% turb vs wave dominance
% orbital velocity
wvnum = 14; % wave number from sine fit
wvlen = 2*pi/wvnum; % wavelength
amp = 0.01; % amplitude from sine fit
T = 0.3; % period from wave gauge spectrum
omega = 2*pi/T; % angular frequency
u_wave = omega*amp/tanh(wvlen*H);


%% stdev of pz vs depth
pz_stdev = cell(size(n));
for i = 1:length(n)
    if nonsphere(i)
        pzmean{i} = condition_vars(abs(smangles{i}(:,2)),smtracks{i}(:,2),Nbins,scaleflag,binedg); 
        pzfluct = abs(smangles{i}(:,2)) - interp1(zprof{i},pzmean{i},smtracks{i}(:,2));
        [pz_stdev{i}, ~] = condition_vars(pzfluct.^2,smtracks{i}(:,2),Nbins,scaleflag,binedg);
    end
end
figure; l=compare_plots(pz_stdev(nonsphere),zprof(nonsphere),mk_col(nonsphere),mk(nonsphere),mk_sz(nonsphere),ls(nonsphere));
xlabel('stdev($p_z$)'); ylabel('$z$ [m]'); legend(l,lstr(nonsphere),'location','se');
goodplot([5 4])


%% concentration as function of wave phase
% Nbins_phs = 6;
% phs_bins_edg = linspace(0,2*pi,Nbins_phs+1);
% Nbins2 = 20;
% 
% lstr_phs = cell(Nbins_phs,1);
% pdf_col = mat2cell(parula(Nbins_phs+1),ones(1,Nbins_phs+1),3);
% pdf_col = pdf_col(1:end-1);
% pdf_mk = {'.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' ...
%     '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.'};
% pdf_mk_sz = 10*ones(Nbins_phs,1);
% pdf_ls = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' ...
%     '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-'};
% 
% % (profiles)
% Cnorm_phs = cell(size(n));
% zprof_phs = cell(size(n));
% for i = 1:length(n)
%     Cnorm_phs{i} = cell(Nbins_phs,1); %zeros(Nbins,Nbins_phs);
%     zprof_phs{i} = cell(Nbins_phs,1);
%     for j = 1:Nbins_phs
%         idx = phs_tracks{i} >= phs_bins_edg(j) & phs_tracks{i} < phs_bins_edg(j+1);
%         lstr_phs{j} = num2str(round(mean([phs_bins_edg(j) phs_bins_edg(j+1)]),2));
%         [~, zprof_phs{i}{j}, C_count_phs] = condition_vars(ones(size(smtracks{i}(idx,2))),smtracks{i}(idx,2),Nbins2,scaleflag,binedg); 
%         C_phs = C_count_phs/(ROI_area*diff(zprof_phs{i}{j}(1:2))*Nt(i));
% 
%         % scaling by uniform concentration
%         if scaleflag
%             normval = length(smtracks{i}(:,2))*real(diff(logspace(log10(binedg(1)),log10(binedg(2)),Nbins2+1)'))/diff(binedg)/(ROI_area*diff(zprof{i}(1:2))*Nt(i));
%         else
%             normval = length(smtracks{i}(:,2))/Nbins2/(ROI_area*diff(zprof{i}(1:2))*Nt(i));
%         end
%         Cnorm_phs{i}{j} = C_phs./normval;
%     end 
% end
% 
% % plot
% figure;
% i = 6;
% l = compare_plots(Cnorm_phs{i}, zprof_phs{i}, pdf_col, pdf_mk, pdf_mk_sz, pdf_ls);
% xlabel('$C$'); ylabel('z [m]'); legend(lstr_phs,'location','se')
% title(lstr{i});
% goodplot([5 4])
% 
% %% depth-integrated concentration vs phase
% Cnorm_phs_int = cell(length(n),1);
% for i = 1:length(n)
%     Cnorm_phs_int{i} = zeros(Nbins_phs,1);
%     for j = 1:Nbins_phs   
%         Cnorm_phs_int{i}(j) = trapz(zprof_phs{i}{j},Cnorm_phs{i}{j});
%     end
% end
% figure; l=compare_plots({phs_bins_edg(1:end-1), phs_bins_edg(1:end-1), phs_bins_edg(1:end-1), phs_bins_edg(1:end-1), phs_bins_edg(1:end-1), phs_bins_edg(1:end-1), phs_bins_edg(1:end-1), phs_bins_edg(1:end-1)}, ...
%     Cnorm_phs_int,mk_col,mk,mk_sz, ls);
% xlabel('$\theta$ [rad]'); ylabel('integrated $C/C_0$'); legend(lstr,'location','northeastoutside')
% goodplot([5 4])
% 
% %% condition C on z-tilde
% binedg(2) = -0.03;
% C_tilde = cell(size(n));
% C_count_tilde = cell(size(n));
% zprof_tilde = cell(size(n));
% zprof_k = cell(size(n)); % z normalized by dom wavenumber
% scaleflag = 0; % code not compatible with scaleflag=1!
% for i = 1:length(n)
%     % mean concentration
%     [~, zprof_tilde{i}, C_count_tilde{i}] = condition_vars(ones(size(smtracks{i}(:,2))),ztilde{i},Nbins2,scaleflag,binedg); 
%     C_tilde{i} = C_count_tilde{i}/(ROI_area*diff(zprof_tilde{i}(1:2))*Nt(i));
% 
%     % scaling by uniform concentration
%     if scaleflag
%         normval = length(smtracks{i}(:,2))*real(diff(logspace(log10(binedg(1)),log10(binedg(2)),Nbins2+1)'))/diff(binedg)/(ROI_area*diff(zprof_tilde{i}(1:2))*Nt(i));
%     else
%         normval = length(smtracks{i}(:,2))/Nbins2/(ROI_area*diff(zprof_tilde{i}(1:2))*Nt(i));
%     end
%     Cnorm_tilde{i} = C_tilde{i}./normval;
%     Cnorm_tilde{i}(Cnorm_tilde{i}==0) = 1e-4;
% end
% 
% % plot concentration
% figure;
% l = compare_plots(Cnorm_tilde, zprof_tilde, mk_col, mk, mk_sz, ls, ...
%     num2cell(nan(size(n))), num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% % set(gca,'yscale','log'); 
% axis([0 4 prof_zlim])
% xlabel('$C/C_0$'); ylabel('$z$ [m]'); legend(l,lstr,'location','se')
% goodplot([6 5])


%% ---------- experimental stuff below ------------ % %
return























%% mean, wave, and turbulent velocity
%% all tracks at once
smtracks_turb = smtracks;
for i = 1:length(n)
%     [E, f] = get_spectrum(smtracks{i}(:,4), fs);
%     figure; loglog(f,E); hold on
%     xlabel('f [Hz]'); ylabel('E_w(f) [m^2/s^2/Hz]')
    
    Yx = fft(smtracks{i}(:,3));
    Yz = fft(smtracks{i}(:,4));
    L = length(Yx);
    f = (1:L)-ceil(L/2);
    f = f*fs/L;
    f = fftshift(f);
    % figure;plot(f,Yx);ylim([-200 200])
    f_co = 2.5;
    Yx(abs(f)<f_co) = 0;%Yx(abs(f)<f_co) - real(Yx(abs(f)<f_co));
    Yz(abs(f)<f_co) = 0;%Yz(abs(f)<f_co) - real(Yz(abs(f)<f_co));
    u_filtered = real(ifft(Yx));
    w_filtered = real(ifft(Yz));
    
%     [E2, f] = get_spectrum(w_filtered, fs);
%     loglog(f,E2); ylim([min(E) max(E)])
%     loglog(f(f>3),1e-5*f(f>3).^(-5/3),'-k')
%     legend('original','filtered','location','sw')
    
%     figure;
%     track_ids = round(linspace(2,length(smtracklength{1}),100));  % 1:30; %
%     c = turbo(100); %turbo(length(track_ids));
%     velmags = u_filtered; %sqrt(u_filtered.^2 + w_filtered.^2);
%     velmags = velmags - min(velmags);
%     velmags = round(velmags*99/max(velmags)) + 1;
%     for j = 1:length(track_ids)
%         idx = smtracks{i}(:,5)==track_ids(j);
%         c_idx = velmags(idx); 
%         scatter(smtracks{i}(idx,1),smtracks{i}(idx,2),4,c(c_idx,:),'filled');
%         hold on
%     end
%     axis equal; axis([-.5 .5 -.45 .05]);
%     xlabel('x [m]'); ylabel('y [m]')
    
%     ufluct_tmp = smtracks{i}(:,4)-mean(smtracks{i}(:,4));
%     t_snip = 1000:1500;
%     figure; plot(t_snip/fs,ufluct_tmp(t_snip),'.-',t_snip/fs,w_filtered(t_snip),'.-'); 
%     xlabel('t [s]'); ylabel('w_p [m/s]'); legend('w_p - \langlew_p\rangle','w_{p,filtered}')
    
    smtracks_turb{i}(:,3) = u_filtered;
    smtracks_turb{i}(:,4) = w_filtered;
end

%% track-by-track
smtracks_turb = smtracks;
for i = 1:length(n)
    for j = 1:max(smtracks{i}(:,5))
        idx = smtracks{i}(:,5) == j;
        
        Yx = fft(smtracks{i}(idx,3));
        Yz = fft(smtracks{i}(idx,4));
        L = length(Yx);
        f = (1:L)-ceil(L/2);
        f = f*fs/L;
        f = fftshift(f);
%         figure;plot(f,Yx);ylim([-200 200])
        f_co = 3;
        Yx(abs(f)<f_co) = 0;%Yx(abs(f)<f_co) - real(Yx(abs(f)<f_co));
        Yz(abs(f)<f_co) = 0;%Yz(abs(f)<f_co) - real(Yz(abs(f)<f_co));
        u_filtered = real(ifft(Yx));
        w_filtered = real(ifft(Yz));
        
        smtracks_turb{i}(idx,3) = u_filtered;
        smtracks_turb{i}(idx,4) = w_filtered;
    end
end
%%
i = 1;
figure;
track_ids = round(linspace(2,length(smtracklength{1}),100));  % 1:30; %
c = turbo(100); %turbo(length(track_ids));
velmags = smtracks_turb{i}(:,3); %sqrt(u_filtered.^2 + w_filtered.^2);
velmags = velmags - min(velmags);
velmags = round(velmags*99/max(velmags)) + 1;
for j = 1:length(track_ids)
    idx = smtracks_turb{i}(:,5)==track_ids(j);
    c_idx = velmags(idx); 
    scatter(smtracks_turb{i}(idx,1),smtracks_turb{i}(idx,2),4,c(c_idx,:),'filled');
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')


%% low-pass filter
kernel = 0;
smtracks_wave = cell(size(n));
smtracklength_wave = cell(size(n));
for i = 1:length(n)
    [smtracks_wave{i}, smtracklength_wave{i}] = smooth_tracks(smtracks{i},kernel,1/fs);
end

% preview low-passed tracks
% track lengths
figure; histogram(smtracklength{1},100)
xlabel('track length [frames]'); ylabel('count')

i = 1;
figure;
track_ids = round(linspace(2,length(smtracklength_wave{i}),100));  % 1:30; %
c = turbo(100); %turbo(length(track_ids));
velmags = smtracks_wave{i}(:,3); %sqrt(smtracks{1}(:,3).^2 + smtracks{1}(:,4).^2);
velmags = velmags - min(velmags);
velmags = round(velmags*99/max(velmags)) + 1;
for j = 1:length(track_ids)
    idx = smtracks_wave{i}(:,5)==track_ids(j);
    c_idx = velmags(idx); % round(smtracklength(track_ids(j))/max(smtracklength(track_ids))*length(track_ids));
    scatter(smtracks_wave{i}(idx,1),smtracks_wave{i}(idx,2),4,c(c_idx,:),'filled');
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')




%% Lagrangian particle vel/accel autocorrelation

% turbulence not yet subtracted
smtracks_turb = smtracks;

% 3 for u, 4 for w, 8 for a_x, 9 for a_z
q_idx = [4];%,8,9];
r_up = cell(length(n),1); 
t_lags = cell(length(n),1);
c = parula(Nbins_wide+1);
c = c(1:end-1,:);
labelstr = {'u_p','w_p','a_{x,p}','a_{z,p}'};

for k = 1:length(n)
    Tmax = max(smtracklength{k});
    r_up{k} = zeros(Nbins_wide,Tmax,length(q_idx));
    
    for i = 1:length(q_idx)
        q1 = q_idx(i); q2 = q1;     
        [r_up{k}(:,:,i), N_pts, t_lags{k}] = lagrang_xcorr(smtracks_turb{k}, q1, q2, z_bins_wide_edg, 0, Tmax, fs, 0);
    
        % plot
        lgnd = cell(Nbins_wide,1); 
        l1 = zeros(Nbins_wide,1);
        figure; 
        Nmin = 100;
        for j = 1:Nbins_wide 
            idx = N_pts(j,:)>Nmin;
            l1(j) = plot(t_lags{k}(idx)/t_plus, r_up{k}(j,idx,i),'.','markersize',10,'color',c(j,:)); hold on
            lgnd{j} = ['$y_0^+ = ' num2str(z_bins_wide_cen(j)) '$'];
        end
        xlabel('$\Delta t^+$'); ylabel(['$\rho_{' labelstr{i} '}(\Delta t)$']); %grid on;
        % set(gca,'YTickLabels',[]);
        line(get(gca,'XLim'),[0 0],'color',[.5 .5 .5])
%         axis([0 80 0 1]); set(gca,'XTick',0:20:100)
        if q1 == 3
            legend(l1,lgnd,'location','northeast'); 
        end
    
        goodplot([4 3.5])

    end
end


%% diffusivity from autocorr
% increase number of Nbins_wide
% check exponential fits on near-surface bins; may need to shorten the
%   fitting range

% 3 for u, 4 for v, 8 for a_x, 9 for a_y, 10 for u_f_at_p, 11 for v_f_at_p, 13 for vtp, 14 for atp, 15 for anp
eps_p_Ta = cell(length(n),1);
for k = 1:length(n)
    q1 = 4; q2 = q1;     
    Tmax = max(smtracklength{k});
    Tint = zeros(Nbins_wide,2);
    
    [r_up, N_pts, t_lags] = lagrang_xcorr(smtracks{k}, q1, q2, z_bins_wide_edg, 0, Tmax, fs, 0);

    Tint = zeros(Nbins_wide,1);
    figure;
    for j = 1:Nbins_wide
        switch q1  % set range for exponential fit (in viscous time units)
            case 4
                t1 = 0; t2 = 4;
        end
        if N_pts(j,t2) > 0
            P = polyfit(t_lags(t_lags<t2 & t_lags>t1 & r_up(j,:)>0), log(r_up(j,t_lags<t2 & t_lags>t1 & r_up(j,:)>0)),1);
            r_up_exp = exp(P(2))*exp(P(1)*t_lags);
            plot(t_lags, r_up_exp,'--','color',c(j,:)); hold on
            plot(t_lags(N_pts(j,:)>Nmin), r_up(j,N_pts(j,:)>Nmin),'.','markersize',10,'color',c(j,:));
            Tint(j) = interp1(r_up_exp,t_lags,1/exp(1));
            xlim([0 10]);
        else
            Tint(j) = nan;
        end
    end

    % vertical diff
    var_wp = zeros(Nbins_wide,1);
    for j = 1:Nbins_wide
        idx_j = smtracks{k}(:,2) >= z_bins_wide_edg(j) & smtracks{k}(:,2) < z_bins_wide_edg(j+1);
        var_wp(j) = var(smtracks{k}(idx_j,q1)-interp1(zprof{k},wmean{k},smtracks{k}(idx_j,2)),'omitnan');
    end
    eps_p_Ta{k} = Tint.*var_wp;
end

figure;
compare_plots(eps_p_Ta,{z_bins_wide_cen,z_bins_wide_cen,z_bins_wide_cen,z_bins_wide_cen,z_bins_wide_cen,z_bins_wide_cen,z_bins_wide_cen,z_bins_wide_cen},...
    mk_col,mk,mk_sz,ls);
if q1 == 4
%             hold on; compare_plots(eps_p_flux,z_prof,'ks','linewidth',1)
%             semilogy(eps_f_RP/nu,yp/del_nu,'k--','linewidth',1)
%             semilogy(eps_f_ML/nu,y/del_nu,'k-','linewidth',1)
    ylabel('$z$ [m]'); ylim(prof_zlim); 
    xlabel('$\varepsilon$'); 
%             legend('$\varepsilon_{p,x}$','$\varepsilon_{p,y}$','$\varepsilon_{p,y}$ (flux)','$\varepsilon_{f}$ (RP)','$\varepsilon_{f}$ (ML)','location','northeastoutside')
    goodplot([5 3.5])
end