% process MP
% process microplastic particle positions/orientations

close all
clear

gdrive_path = 'C:\Users\ljbak\My Drive\';  % 'H:\My Drive\';  %   'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set
n = [15,13,1:6]; %  runs to include
% n = [16,14,7:12];

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\run_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fs = run_params.imagingFreq_Hz(1);

nonsphere = zeros(size(n));
for i = 1:length(n)
    nonsphere(i) = strncmp(run_params.ParticleType{n(i)},'d',1) || strncmp(run_params.ParticleType{n(i)},'r',1);
end

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
    end
end

% % flip streamwise coord so that flow moves in +x direction (flip y to preserve righthandedness)
% smtracks(:,[1 3]) = -smtracks(:,[1 3 8]);
% smangles(:,[2]) = -smangles(:,[2]);             % (?)
% smangles_cont(:,[2]) = -smangles_cont(:,[2]);             % (?)

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
ebar_gray = [.4 .4 .4]; 
ebar_red = [1 .3 .3]; 
ebar_blue = [.3 .3 1]; 
ebar_green = [0.4 0.7 0.4];
sym_green = [0 .5 0];

mk_col = {'k' 'k' 'k' 'r' 'r' 'b' 'b' 'b'};
mk = {'.' '.' '.' '+' '+' 'o' 'o' 'o'};
mk_sz = [6 9 12 5 10 3 5 7];
ls = {'-' '-' '-' '-' '-' '-' '-' '-'};
eb_col = {ebar_gray ebar_gray ebar_gray ebar_red ebar_red ebar_blue ebar_blue ebar_blue};
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
case_idx = 2+2;
load(sprintf('centers_run%02d.mat',n(case_idx)),'z_freesurf_inst_rect');
%%
figure; 
set(gcf,'position',[369.8000  172.2000  817.6000  512.8000]); %[0.0010    0.0410    1.5360    0.7488]*1000)
track_ids = round(linspace(2,length(smtracklength{1}),100));  % 1:30; %
c = turbo(100); %turbo(length(track_ids));
velmags = smtracks{case_idx}(:,4); %sqrt(smtracks{1}(:,3).^2 + smtracks{1}(:,4).^2); % 
velmags = velmags - min(velmags);
velmags = round(velmags*99/max(velmags)) + 1;
stay_time = 60; % number of frames to keep on figure
for i = 1%:500%max(smtracks{1}(:,7))
    idx1 = smtracks{case_idx}(:,7)==i; 
    idx2 = smtracks{case_idx}(:,7)<i & smtracks{case_idx}(:,7)>(i-stay_time);
    c_idx1 = velmags(idx1); %round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));  %   
    c_idx2 = velmags(idx2); 
    scatter(smtracks{case_idx}(idx1,1),smtracks{case_idx}(idx1,2),10,c(c_idx1,:),'filled'); hold on
    scatter(smtracks{case_idx}(idx2,1),smtracks{case_idx}(idx2,2),2,c(c_idx2,:),'filled'); 
    plot(-z_freesurf_inst_rect{i}(:,1),z_freesurf_inst_rect{i}(:,2),'k-','linewidth',0.5); 
    hold off
    axis equal; axis([-.5 .5 -.45 .05]);
    xlabel('x [m]'); ylabel('y [m]')
    goodplot([7 5])

    pause(1/10)
%     fig_to_gif('figs/tracks-r10-16.gif',0.1,'-r600')
end


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
w_C = cell(size(n));
zprof = cell(size(n));
zprof_k = cell(size(n)); % z normalized by dom wavenumber
scaleflag = 0; % code not compatible with scaleflag=1!
for i = 1:length(n)
    % mean concentration
    [~, zprof{i}, C{i}] = condition_vars(ones(size(smtracks{i}(:,2))),smtracks{i}(:,2),Nbins,scaleflag,binedg); 
    C{i} = C{i}/(ROI_area*diff(zprof{i}(1:2))*Nt(i));
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
l = compare_plots(Cnorm, zprof, mk_col, mk, mk_sz, ls, ...
    w_C, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% set(gca,'yscale','log'); 
axis([0 4 prof_zlim])
xlabel('$C/C_0$'); ylabel('$z$ [m]'); legend(l,lstr,'location','se')
goodplot([6 5])

%% theoretical profile
dz = 0.07; % 0.055; % 
Lm_fit = zeros(size(n));
C0 = zeros(size(n));
for i = 1:length(n)
    P = polyfit((zprof{i}+dz),log(Cnorm{i}),1);
    Lm_fit(i) = 1/P(1);
    C0(i) = exp(P(2));
end

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
    z_Lm{i} = (zprof{i}+dz)/Lm_fit(i); % Lm_theor(i); %
end
figure; l = compare_plots(C_C0, z_Lm, mk_col, mk, mk_sz, ls, ...
    w_C_C0, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
hold on; l2 = plot(exp(-8:.1:0),-8:.1:0,'k--','linewidth',2);
ylim([-3 0])
xlabel('$C/C_0$'); ylabel('$(z-\Delta z)/L_m$'); legend([l;l2],[lstr,{'Fit'}],'location','se')
goodplot([6 5])

fprintf('\nLm = %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f m\n',Lm_fit)
fprintf('C0 = %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f\n',C0)

% % check if developed 
% x_us = .3; x_ds = -x_us; us_idx = smtracks{1}(:,1) > x_us; ds_idx = smtracks{1}(:,1) < x_ds;
% [~, ~, C_us] = condition_vars(ones(size(smtracks{1}(us_idx,2))),smtracks{1}(us_idx,2),Nbins,scaleflag,binedg);
% [~, ~, C_ds] = condition_vars(ones(size(smtracks{1}(ds_idx,2))),smtracks{1}(ds_idx,2),Nbins,scaleflag,binedg);
% figure; plot(C_us, zprof{1}, 'b.', C_ds, zprof{1}, 'r.'); xlabel('C'); ylabel('z'); legend('upstream','downstream','location','se')


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
A0 = 1.5*[1 2 3]*1e-2*0.4*5e-2;

% % plot flux balance components
% figure; l = compare_plots(WC,zprof_k,mk_col,mk,mk_sz,ls);
% xlabel('$W_b\langle C\rangle$'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([5 4])
% figure; l = compare_plots(dCdz,zprof_k,mk_col,mk,mk_sz,ls);
% xlabel('$d\langle C\rangle dz$'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([5 4])
% figure; l = compare_plots(fluxz,zprof_k,mk_col,mk,mk_sz,ls);
% xlabel('$\Phi$ [m$^{-2}$s$^{-1}$]'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([5 4])

% plot diffusivity
figure; l = compare_plots(eps_p_flux,zprof_k,mk_col,mk,mk_sz,ls);
hold on; line([A0(1) A0(1)],binedg*kdom,'color','k','linestyle','--'); 
line([A0(2) A0(2)],binedg*kdom,'color','k','linestyle','--'); 
line([A0(3) A0(3)],binedg*kdom,'color','k','linestyle','--'); 
xlim([-.01 .02]); xlabel('$A(z)$ [m$^{2}$s$^{-1}$]'); ylabel('$z k_{dom}$'); legend(l,lstr,'location','se'); goodplot([4 4])



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

% plot
figure;
compare_plots(delta_t_rng, delta_t_pdf, mk_col, mk, mk_sz, ls);
xlabel('$\Delta t$ [s]'); ylabel('PDF'); set(gca,'yscale','log')
goodplot([5 4])

figure;
compare_plots(delta_zmax_rng, delta_zmax_pdf, mk_col, mk, mk_sz, ls);
xlabel('$\Delta z_{max}$ [m]'); ylabel('PDF'); set(gca,'yscale','log')
goodplot([5 4])



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
pdf_col = mat2cell(parula(Nbins_wide),ones(1,Nbins_wide),3);
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

% plot velocity profiles
figure;
subplot(121); 
l = compare_plots(umean, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_umean, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle u_p\rangle$ [m/s]'); ylabel('$zk_{dom}$')
legend(l,lstr,'location','sw')

subplot(122); 
compare_plots(wmean, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_wmean, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle w_p\rangle$ [m/s]'); 
goodplot([5 4])

figure;
subplot(131); 
l = compare_plots(uu, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_uu, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle u_p''u_p''\rangle$ [m$^2$/s$^2$]'); ylabel('$zk_{dom}$')
legend(l,lstr,'location','se')
xlim([0 10]*1e-3)

subplot(132); 
compare_plots(ww, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_ww, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle w_p''w_p''\rangle$ [m$^2$/s$^2$]'); 
xlim([0 10]*1e-3)

subplot(133); 
compare_plots(uw, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_uw, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle u_p''w_p''\rangle$ [m$^2$/s$^2$]'); 
xlim([-5 5]*1e-3)
goodplot([6 4])

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
    
        for j = 1:Nbins_wide
            idx = smtracks{i}(:,2) >= z_bins_wide_edg(j) & smtracks{i}(:,2) < z_bins_wide_edg(j+1);
            [px_pdf{i}(:,j), px_rng{i}(:,j)] = pdf_var(abs(smangles{i}(idx,1)),Npdf,0,[0 1]);
            [py_pdf{i}(:,j), py_rng{i}(:,j)] = pdf_var(abs(smangles{i}(idx,3)),Npdf,0,[0 1]);
            [pz_pdf{i}(:,j), pz_rng{i}(:,j)] = pdf_var(abs(smangles{i}(idx,2)),Npdf,0,[0 1]);
            [pzs_pdf{i}(:,j), pzs_rng{i}(:,j)] = pdf_var(smangles{i}(idx,2),Npdfs,0,[-1 1]);
        end
    end
end

% plot orientation pdfs
pdf_col = mat2cell(parula(Nbins_wide),ones(1,Nbins_wide),3);
pdf_mk = {'.' '.' '.' '.' '.' '.'};
pdf_mk_sz = 10*ones(Nbins_wide,1);
pdf_ls = {'-' '-' '-' '-' '-' '-'};

for i = 1:length(n)
    if nonsphere(i)
        figure;
        subplot(141);
        l = compare_plots(mat2cell(px_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(px_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
        xlabel('$|p_x|$'); ylabel('PDF')
        subplot(142);
        compare_plots(mat2cell(py_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(py_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
        xlabel('$|p_y|$')
        subplot(143);
        compare_plots(mat2cell(pz_rng{i},Npdf,ones(1,Nbins_wide)), ...
            mat2cell(pz_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
        xlabel('$|p_z|$')
        subplot(144);
        compare_plots(mat2cell(pzs_rng{i},Npdfs,ones(1,Nbins_wide)), ...
            mat2cell(pzs_pdf{i},Npdfs,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
        xlabel('$p_z$')
        sgtitle(run_params.ParticleType{n(i)})
        goodplot([6 4])
    end
end

%% plot pz pdfs only
pdf_col = mat2cell(parula(Nbins_wide),ones(1,Nbins_wide),3);
pdf_mk = {'.' '.' '.' '.' '.' '.'};
pdf_mk_sz = 10*ones(Nbins_wide,1);
pdf_ls = {'-' '-' '-' '-' '-' '-'};

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


%% depth & orientation = irradiation

% radiation as a function of depth 
Lz = .2; % light decay lengthscale [m]  (ARBITRARY)
rad_level = @(z) exp(z/Lz);   % normalized by surface intensity
% rad_level = @(z) 1;

Anorm = cell(size(n));      % particle planar area normal to vertical, not normalized [m^2]
Anorm_max = cell(size(n));  % particle planar area normal to vertical normalized by MAX AREA
Anorm_sph = cell(size(n));  % particle planar area normal to vertical normalized by EQUIVALENT SPHERE AREA
irrad = cell(size(n));      % irradiation [W]
irrad_max = cell(size(n));  % irradiation on particle normalized by its max planar area and surface light level
irrad_sph = cell(size(n));  % irradiation on particle normalized by its equiv sphere planar area and surface light level
irrad_mean = zeros(size(n));
irrad_max_mean = zeros(size(n));
irrad_sph_mean = zeros(size(n));

for i = 1:length(n)
    Anorm{i} = ones(length(smtracks{i}),1)*run_params.Dp_m(n(i))^2/4*pi;
    Anorm_max{i} = ones(length(smtracks{i}),1);
    Anorm_sph{i} = ones(length(smtracks{i}),1);
    if nonsphere(i) 
        if strncmp(run_params.ParticleType{n(i)},'d',1)
            Anorm{i} = abs(smangles{i}(:,2))*run_params.Dp_m(n(i))^2/4*pi;
            Anorm_max{i} = abs(smangles{i}(:,2));  
            Anorm_sph{i} = Anorm_max{i}/(3/2*AR(i))^(2/3);
        else
            Anorm{i} = abs(smangles{i}(:,2))*run_params.Dp_m(n(i))*1.75e-3;
            Anorm_max{i} = sqrt(1 - (smangles{i}(:,2).^2));
            Anorm_sph{i} = Anorm_max{i}/(pi*(9/16*1/AR(i))^(1/3));
        end
    end
    irrad{i} = Anorm{i}.*rad_level(smtracks{i}(:,2));
    irrad_max{i} = Anorm_max{i}.*rad_level(smtracks{i}(:,2));
    irrad_sph{i} = Anorm_sph{i}.*rad_level(smtracks{i}(:,2));
    irrad_mean(i) = mean(irrad{i},'omitnan');
    irrad_max_mean(i) = mean(irrad_max{i},'omitnan');
    irrad_sph_mean(i) = mean(irrad_sph{i},'omitnan');
end

fprintf(['irradiation normalized by surface intensity: [m^-2]\n' ...
    '\t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \n'], lstr{:});
disp(irrad_mean)
fprintf(['irradiation normalized by max planar area and surface intensity:\n' ...
    '\t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \n'], lstr{:});
disp(irrad_max_mean)
fprintf(['irradiation normalized by planar area of equiv sphere and surface intensity:\n' ...
    '\t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \t\t%s \n'], lstr{:})
disp(irrad_sph_mean)

% plot
figure;
compare_plots(num2cell(run_params.riseVel_m_s(n)), num2cell(irrad_max_mean),mk_col,mk,mk_sz,ls);
xlabel('W_b [m/s]'); ylabel('I_{norm,max}')

figure;
compare_plots(num2cell(run_params.riseVel_m_s(n)), num2cell(irrad_sph_mean),mk_col,mk,mk_sz,ls);
xlabel('W_b [m/s]'); ylabel('I_{norm,sph}')

figure;
compare_plots(num2cell(Lm_fit), num2cell(irrad_max_mean),mk_col,mk,mk_sz,ls);
xlabel('L_m [m]'); ylabel('I_{norm,max}')

figure;
compare_plots(num2cell(Lm_fit), num2cell(irrad_sph_mean),mk_col,mk,mk_sz,ls);
xlabel('L_m [m]'); ylabel('I_{norm,sph}')


% % ---------- experimental stuff below ------------ % %
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
c = parula(Nbins_wide);
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