% process MP
% process microplastic particle positions/orientations

close all
clear

gdrive_path = 'C:\Users\ljbak\My Drive\';  %  'H:\My Drive\';  %  'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220315';  % expt set
n = 1:6; %  7:12; %  runs to include

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

% compute error bars (90% CI)
T_indep = 5;  % number of frames per independent realization
error_mean = @(q_var,q_N) 1.645*sqrt(q_var./(q_N/T_indep));
error_var = @(q_var,q_N) q_var.*(1-(q_N/T_indep-1)./chi2inv(0.05,q_N/T_indep-1));

% binning vars
Nbins = 12;   % number of bins in profiles
binedg = [-.45 -.05];
prof_ylim = [-.5 0];

Nbins_wide = 4;
y_bins_wide_edg = real(logspace(log10(binedg(1)),log10(binedg(2)),Nbins_wide+1)); % wide bin edges for autocorrelations and PDFs
y_bins_wide_cen = mean([y_bins_wide_edg(2:end); y_bins_wide_edg(1:end-1)],1);

Nt = zeros(size(n));
for i = 1:length(n); Nt(i) = max(smtracks{i}(:,7)); end

% plotting vars
ebar_gray = [.4 .4 .4]; 
ebar_red = [1 .3 .3]; 
ebar_blue = [.3 .3 1]; 
ebar_green = [0.4 0.7 0.4];
sym_green = [0 .5 0];

mk_col = {'k' 'r' 'r' 'b' 'b' 'b'};
mk = {'.' '+' '+' 'o' 'o' 'o'};
mk_sz = [10 5 10 3 5 7];
ls = {'-' '-' '-' '-' '-' '-'};
eb_col = {ebar_gray ebar_red ebar_red ebar_blue ebar_blue ebar_blue};
lstr = {'n','r10','r20','d5','d7','d10'};

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
    if strncmp(run_params.ParticleType{i},'d',1)
        AR(i) = 1.6e-3/run_params.Dp_m(i);
    elseif strncmp(run_params.ParticleType{i},'r',1)
        AR(i) = run_params.Dp_m(i)/1.3e-3;
    end
end


%% preview tracks
% track lengths
figure; histogram(smtracklength{1},100)
xlabel('track length [frames]'); ylabel('count')

figure;
track_ids = round(linspace(2,length(smtracklength{1}),100));  % 1:30; %
c = jet(100); %jet(length(track_ids));
velmags = smtracks{1}(:,3); %sqrt(smtracks{1}(:,3).^2 + smtracks{1}(:,4).^2);
velmags = velmags - min(velmags);
velmags = round(velmags*99/max(velmags)) + 1;
for i = 1:length(track_ids)
    idx = smtracks{1}(:,5)==track_ids(i);
    c_idx = velmags(idx); % round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));
    scatter(smtracks{1}(idx,1),smtracks{1}(idx,2),4,c(c_idx,:),'filled');
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')


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
    w_C{i} = w_C{i}./normval;
end

% plot concentration
figure;
l = compare_plots(Cnorm, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_C, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
% set(gca,'yscale','log'); 
axis([0 5 prof_ylim*kdom])
xlabel('$C/C_0$'); ylabel('$zk_w$'); legend(l,lstr,'location','se')
goodplot([5 4])

%% theoretical profile
Lm_fit = zeros(size(n));
C0 = zeros(size(n));
for i = 1:length(n)
    P = polyfit(zprof{i},log(Cnorm{i}),1);
    Lm_fit(i) = 1/P(1);
    C0(i) = exp(P(2));
end

A0 = 1.5*ut*0.4*Hs;
Lm_theor = A0./run_params.riseVel_m_s(n)';

% i = 1;
% z_exp = linspace(zprof{i}(1),zprof{i}(end),50);
% C_exp = C0(i)*exp(z_exp/Lm_fit(i));
% hold on; plot(C_exp,z_exp*kdom,'k--','linewidth',1)


%% collapse conc profile
C_C0 = Cnorm;
w_C_C0 = w_C;
z_Lm = zprof;
for i = 1:length(n)
    C_C0{i} = Cnorm{i}/C0(i);
    w_C_C0{i} = w_C{i}/C0(i);
    z_Lm{i} = zprof{i}/Lm_theor(i); %Lm_fit(i);
end
figure; l = compare_plots(C_C0, z_Lm, mk_col, mk, mk_sz, ls, ...
    w_C_C0, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
hold on; plot(exp(-8:.1:0),-8:.1:0,'k--')
xlabel('$C/C_0$'); ylabel('$z/L_m$'); legend(l,lstr,'location','se')
goodplot([5 4])

fprintf('\nLm = %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f m\n',Lm_fit)
fprintf('C0 = %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f\n',C0)




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
dCdz = cell(size(n));
WC = cell(size(n));
for i = 1:length(n)
    dCdz{i} = gradient(C{i})./diff(zprof{i}(1:2));
    WC{i} = run_params.riseVel_m_s(n(i))*C{i};
end

% solve for diffusivity
W_rise = run_params.riseVel_m_s;
Az = cell(size(n));
for i = 1:length(n)
    Az{i} = (WC{i} - fluxz{i})./dCdz{i};
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
figure; l = compare_plots(Az,zprof_k,mk_col,mk,mk_sz,ls);
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
        idx = smtracks{i}(:,2) >= y_bins_wide_edg(j) & smtracks{i}(:,2) < y_bins_wide_edg(j+1);
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

subplot(132); 
compare_plots(ww, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_ww, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle w_p''w_p''\rangle$ [m$^2$/s$^2$]'); 

subplot(133); 
compare_plots(uw, zprof_k, mk_col, mk, mk_sz, ls, ...
    w_uw, num2cell(nan(size(n))), eb_col, num2cell(nan(size(n))));
xlabel('$\langle u_p''w_p''\rangle$ [m$^2$/s$^2$]'); 
goodplot([6 4])


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
            idx = smtracks{i}(:,2) >= y_bins_wide_edg(j) & smtracks{i}(:,2) < y_bins_wide_edg(j+1);
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

% for i = 1:length(n)
%     if nonsphere(i)
%         figure;
%         subplot(141);
%         l = compare_plots(mat2cell(px_rng{i},Npdf,ones(1,Nbins_wide)), ...
%             mat2cell(px_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%         xlabel('$|p_x|$'); ylabel('PDF')
%         subplot(142);
%         compare_plots(mat2cell(py_rng{i},Npdf,ones(1,Nbins_wide)), ...
%             mat2cell(py_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%         xlabel('$|p_y|$')
%         subplot(143);
%         compare_plots(mat2cell(pz_rng{i},Npdf,ones(1,Nbins_wide)), ...
%             mat2cell(pz_pdf{i},Npdf,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%         xlabel('$|p_z|$')
%         subplot(144);
%         compare_plots(mat2cell(pzs_rng{i},Npdfs,ones(1,Nbins_wide)), ...
%             mat2cell(pzs_pdf{i},Npdfs,ones(1,Nbins_wide)), pdf_col,pdf_mk,pdf_mk_sz,pdf_ls);
%         xlabel('$p_z$')
%         sgtitle(run_params.ParticleType{i})
%         goodplot([6 4])
%     end
% end

%% depth & orientation = irradiation

% radiation as a function of depth 
Lz = -.1; % light decay lengthscale [m]
rad_level = @(z) exp(z/Lz);   % (CHECK THIS)

Anorm_max = cell(size(n));  % particle planar area normal to vertical normalized by MAX AREA
Anorm_sph = cell(size(n));  % particle planar area normal to vertical normalized by EQUIVALENT SPHERE AREA
irrad = cell(size(n));  % irradiation on particle normalized by its nominal planar area and surface light level

for i = 1:length(n)
    Anorm_max{i} = ones(length(smtracks{i}),1);
    Anorm_sph{i} = ones(length(smtracks{i}),1);
    if nonsphere(i) 
        if strncmp(run_params.ParticleType,'d',1)
            Anorm_max{i} = abs(smangles{i}(:,2));  
            Anorm_sph{i} = Anorm_max{i}/(3/2*AR(i))^(2/3);
        else
            Anorm_max{i} = sqrt(1 - (smangles{i}(:,2).^2));
            Anorm_sph{i} = Anorm_max{i}/(pi*(9/16*1/AR(i))^(1/3));
        end
    end
    irrad_max{i} = Anorm_max{i}.*rad_level(smtracks{i}(:,2));
    irrad_sph{i} = Anorm_sph{i}.*rad_level(smtracks{i}(:,2));
    disp(mean(irrad_max{i},'omitnan'))
%     disp(mean(irrad_sph{i},'omitnan'))
end




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
%     c = jet(100); %jet(length(track_ids));
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
    for j = 1:length(smtracklength)
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







%% Lagrangian particle vel/accel autocorrelation

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
        [r_up{k}(:,:,i), N_pts, t_lags{k}] = lagrang_xcorr(smtracks_turb{k}, q1, q2, y_bins_wide_edg, 0, Tmax, fs, 0);
    
        % plot
        lgnd = cell(Nbins_wide,1); 
        l1 = zeros(Nbins_wide,1);
        figure; 
        Nmin = 100;
        for j = 1:Nbins_wide 
            idx = N_pts(j,:)>Nmin;
            l1(j) = plot(t_lags{k}(idx)/t_plus, r_up{k}(j,idx,i),'.','markersize',10,'color',c(j,:)); hold on
            lgnd{j} = ['$y_0^+ = ' num2str(mean(y_bins_wide_edg(j:j+1)),2) '$'];
        end
        xlabel('$\Delta t^+$'); ylabel(['$\rho_{' labelstr{i} '}(\Delta t)$']); %grid on;
        % set(gca,'YTickLabels',[]);
        line(get(gca,'XLim'),[0 0],'color',[.5 .5 .5])
%         axis([0 80 0 1]); set(gca,'XTick',0:20:100)
        if q1 == 3
            legend(l1,lgnd,'location','northeast'); 
        end
    
        goodplot([4 3.5])
        if q1 == 8
            goodplot([5 3.5]);
%             axis([0 50 -.25 1])
            legend(l1,lgnd,'location','northeast'); 
        end
    end
end

