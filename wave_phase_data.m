% wave phase data
load('run_1_dewarped.mat');

%%

tt = 1/30 * (1:(size(run_dewarp,2)/4));
figure; plot(tt,run_dewarp(100,1:end/4));
xlabel('t [s]'); ylabel('z [m]')
xlim([0 5])

zz = 1:size(run_dewarp,1);
figure; plot(zz,run_dewarp(:,100));
xlabel('x [px]'); ylabel('z [m]')
xlim([0 5])

% % Hilbert transform
% x = hilbert(run_dewarp(1000,:));
% phs = atan2(imag(x),real(x));

