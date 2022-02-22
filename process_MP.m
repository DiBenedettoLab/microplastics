% process MP
% process microplastic particle positions/orientations

close all
clear

% file paths of runs to include (for real expt: load from spreadsheet)
fpaths = {'wasirf particle test 6 - nurdles/', ...
    'wasirf particle test 7 - small rods/', ...
    'wasirf particle test 8 - small disks/'};

s = load([fpaths{1} 'centers.mat']);
r = load([fpaths{2} 'centers.mat']);
d = load([fpaths{3} 'centers.mat']);

pxtom = 5e-3/(2048*.29/21.33);

s.zp = [];
for i = 1:length(s.centers)
    s.zp = [s.zp; s.centers{i}(:,2)];
end
r.zp = [];
for i = 1:length(r.centers)
    r.zp = [r.zp; r.centers{i}(:,2)];
end
d.zp = [];
for i = 1:length(d.centers)
    d.zp = [d.zp; d.centers{i}(:,2)];
end

s.N_bins = 18;
r.N_bins = 15;
d.N_bins = 12;
s.lims = linspace(min(s.zp),max(s.zp),s.N_bins+1);
r.lims = linspace(min(r.zp),max(r.zp),r.N_bins+1);
d.lims = linspace(min(d.zp),max(d.zp),d.N_bins+1);

ebar_mean = @(q_rms,q_N) 1.645*q_rms./sqrt(q_N);

ebar_gray = [.4 .4 .4]; 
ebar_red = [1 .3 .3]; 
ebar_blue = [.3 .3 1]; 

[~, s.z_prof, s.C] = condition_vars(ones(size(s.zp)),s.zp,s.N_bins,0);
[~, r.z_prof, r.C] = condition_vars(ones(size(r.zp)),r.zp,r.N_bins,0);
[~, d.z_prof, d.C] = condition_vars(ones(size(d.zp)),d.zp,d.N_bins,0);

% uncertainty
C = zeros(length(s.centers),s.N_bins);
for i = 1:length(s.centers)
    zp = s.centers{i}(:,2);
    if ~isempty(zp)
        [~,~,C(i,:)] = condition_vars(ones(size(zp)),zp,s.N_bins,0);
    end
end
s.w_C = ebar_mean(rms(C,1)'*length(s.centers),s.C/5);

C = zeros(length(r.centers),r.N_bins);
for i = 1:length(r.centers)
    zp = r.centers{i}(:,2);
    if ~isempty(zp)
        [~,~,C(i,:)] = condition_vars(ones(size(zp)),zp,r.N_bins,0);
    end
end
r.w_C = ebar_mean(rms(C,1)'*length(r.centers),r.C/5);

C = zeros(length(d.centers),d.N_bins);
for i = 1:length(d.centers)
    zp = d.centers{i}(:,2);
    if ~isempty(zp)
        [~,~,C(i,:)] = condition_vars(ones(size(zp)),zp,d.N_bins,0);
    end
end
d.w_C = ebar_mean(rms(C,1)'*length(d.centers),d.C/5);

% scaling
k_dom = 1/.6;
surf_lvl = (21.33-2)/21.33*2048;
s.z_prof = s.z_prof - surf_lvl;
s.z_prof = s.z_prof*pxtom; 
r.z_prof = r.z_prof - surf_lvl;
r.z_prof = r.z_prof*pxtom; 
d.z_prof = d.z_prof - surf_lvl;
d.z_prof = d.z_prof*pxtom; 

s.C = s.C/(length(s.zp)/s.N_bins); s.w_C = s.w_C/(length(s.zp)/s.N_bins);
r.C = r.C/(length(r.zp)/r.N_bins); r.w_C = r.w_C/(length(r.zp)/r.N_bins);
d.C = d.C/(length(d.zp)/d.N_bins); d.w_C = d.w_C/(length(d.zp)/d.N_bins);

figure;
l = compare_plots({s.C, r.C, d.C}, {s.z_prof*k_dom, r.z_prof*k_dom, d.z_prof*k_dom}, {'k' 'r' 'b'}, {'none' 'none' 'none'}, ...
    [6 6 6], {'-' '-' '-'}, {s.w_C, r.w_C, d.w_C}, {nan nan nan}, {ebar_gray ebar_red ebar_blue}, {nan nan nan});
% l = compare_plots({s.C}, {s.z_prof*k_dom}, {'k'}, {'none'}, ...
%     [6], {'-'}, {s.w_C}, {nan}, {ebar_gray}, {nan});

% theoretical profile
Lz = .2;
z_exp = linspace(-.6,0,50);
C_exp = s.C(end)*exp(z_exp/Lz);
hold on; plot(C_exp,z_exp*k_dom,'k--','linewidth',1)

axis([0 2.5 -.6 0])
xlabel('$C/C_0$'); ylabel('$zk_w$'); 
% legend('Nurdles','location','se')
legend(l,'Nurdles','Rods','Disks','location','se')
goodplot([5 4])
