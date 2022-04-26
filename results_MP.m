% process MP
% process microplastic particle positions/orientations

close all
clear

gdrive_path = 'C:\Users\ljbak\My Drive\';  %  'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220315';  % expt set
n = 1:6; % 7:12; % runs to include

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
for i = 1:(length(n))
    sm_tmp = load(sprintf('smtracks_run%02d.mat',n(i)));
    smtracks{i} = sm_tmp.smtracks;
    smtracklength{i} = sm_tmp.smtracklength;
end

% compute error bars (90% CI)
ebar_mean = @(q_rms,q_N) 1.645*q_rms./sqrt(q_N);

ebar_gray = [.4 .4 .4]; 
ebar_red = [1 .3 .3]; 
ebar_blue = [.3 .3 1]; 
ebar_green = [0.4 0.7 0.4];
sym_green = [0 .5 0];

Nbins = 12;   % number of bins in profiles
Nbins_wide = 5;
binedg = [-.45 -.05];
prof_ylim = [-.8 0];

y_bins_ac = real(logspace(log10(binedg(1)),log10(binedg(2)),Nbins_wide+1)); % wide bins for autocorrelations

T_indep = 5;  % number of frames per independent realization
kdom = 1/.6;  % scaling (dominant wavenumber)
t_plus = 1;



%% concentration

C = cell(size(n));
w_C = cell(size(n));
zprof = cell(size(n));
zprof_k = cell(size(n)); % z normalized by dom wavenumber
scaleflag = 1;
for i = 1:length(n)
    % mean concentration
    [~, zprof{i}, C{i}] = condition_vars(ones(size(smtracks{i}(:,2))),smtracks{i}(:,2),Nbins,scaleflag,binedg); 
    zprof_k{i} = zprof{i}*kdom;
    
    % uncertainty
%     C_tr = zeros(max(smtracks{i}(:,7)),Nbins);
%     for j = 1:max(smtracks{i}(:,7))
%         idx = smtracks{i}(:,7) == j;
%         zp_j = smtracks{i}(idx,2);
%         if ~isempty(zp_j)
%             [~,~,C_tr(j,:)] = condition_vars(ones(size(zp_j)),zp_j,Nbins,0);
%         end
%     end
    w_C{i} = nan; %ebar_mean(rms(C_tr,1)'*length(smtracks{i}(:,2)),C{i}/T_indep);

    % scaling
    if scaleflag
        normval = length(smtracks{i}(:,2))*real(diff(logspace(log10(binedg(1)),log10(binedg(2)),Nbins+1)'))/diff(binedg);
    else
        normval = length(smtracks{i}(:,2))/Nbins;
    end
    C{i} = C{i}./normval;
    w_C{i} = w_C{i}./normval;
end

% plot concentration
figure;
l = compare_plots(C, zprof_k, {'k' 'r' 'r' 'b' 'b' 'b'}, {'.' '+' '+' 'o' 'o' 'o'}, [10 5 10 5 7.5 10], {'-' '-' '-' '-' '-' '-'}, ...
    w_C, num2cell(nan(size(n))), {ebar_gray ebar_red ebar_red ebar_blue ebar_blue ebar_blue}, num2cell(nan(size(n))));

% theoretical profile
Lz = .11;
z_exp = linspace(prof_ylim(1),prof_ylim(2),50);
C_exp = C{1}(end)*exp(z_exp/Lz);
hold on; plot(C_exp,z_exp*kdom,'k--','linewidth',1)

% set(gca,'yscale','log'); 
axis([0 7 prof_ylim*kdom])
xlabel('$C/C_0$'); ylabel('$zk_w$'); 
% legend('Nurdles','location','se')
legend(l,'n','r10','r20','d5','d7','d10','location','se')
goodplot([5 4])


%% mean, wave, and turbulent velocity

% for i = 1:length(smtracks)
idx = smtracks{1}(:,5)==1;
    [E_u, f] = get_spectrum(smtracks{1}(:,3), fs);
    figure; loglog(f,E_u);


%% Lagrangian particle vel/accel autocorrelation

% 3 for u, 4 for w, 8 for a_x, 9 for a_z
q_idx = [3,4];%,8,9];
r_up = cell(length(n),1); 
t_lags = cell(length(n),1);
c = parula(Nbins_wide);
labelstr = {'u_p','w_p','a_{x,p}','a_{z,p}'};

for k = 1:length(n)
    Tmax = max(smtracklength{k});
    r_up{k} = zeros(Nbins_wide,Tmax,length(q_idx));
    
    for i = 1:length(q_idx)
        q1 = q_idx(i); q2 = q1;     
        [r_up{k}(:,:,i), N_pts, t_lags{k}] = lagrang_xcorr(smtracks{k}, q1, q2, y_bins_ac, 0, Tmax, fs, 0);
    
        % plot
        lgnd = cell(Nbins_wide,1); 
        l1 = zeros(Nbins_wide,1);
        figure; 
        Nmin = 100;
        t_lags_long = (0:(2*length(t_lags{k})))/fs;
        for j = 1:Nbins_wide 
            idx = N_pts(j,:)>Nmin;
            l1(j) = plot(t_lags{k}(idx)/t_plus, r_up{k}(j,idx,i),'.','markersize',10,'color',c(j,:)); hold on
            lgnd{j} = ['$y_0^+ = ' num2str(mean(y_bins_ac(j:j+1)),2) '$'];
        end
        xlabel('$\Delta t^+$'); ylabel(['$\rho_{' labelstr{i} '}(\Delta t)$']); %grid on;
        % set(gca,'YTickLabels',[]);
        line(get(gca,'XLim'),[0 0],'color',[.5 .5 .5])
%         axis([0 80 0 1]); set(gca,'XTick',0:20:100)
        if q1 == 3
            legend(l1,lgnd,'location','southwest'); 
        end
    
        goodplot([4 3.5])
        if q1 == 8
            goodplot([5 3.5]);
%             axis([0 50 -.25 1])
            legend(l1,lgnd,'location','northeast'); 
        end
    end
end

