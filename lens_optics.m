% required ray angles
theta1 = 62.5;
theta2 = -45;

% estimate S1 (distance from LED to lens)
S1 = 10e-3; % 20e-3; % 

% calculate best lens focal length using optical invariant and lens-makers eqn
M = tand(theta1)/tand(theta2);
S2 = M*S1;
f = 1/(1/S1 + 1/S2);

fprintf('M = %1.2f, S2 = %1.3f m, f = %1.3f m\n', M, S2, f)

%% choose available focal length
f = 25e-3; %60e-3; % 

% calculate S1 and required lens diameter d
S1 = f*(M+1)/M;
d = 2*S1*tand(theta1);

% available lens diameter and % light transmitted
d_avail = 25.4e-3;
lighttrans = d_avail^2/d^2;

fprintf('f = %1.3f m, S1 = %1.3f m, d = %1.3f m, light transmitted = %2.1f%%\n', f, S1, d, lighttrans*100)

%% vary S1 and observe response on theta2
S1 = linspace(1,50,50)*1e-3;
S2 = 1./(1/f - 1./S1);

% calculate theta2
theta2 = atand(tand(theta1)./(S2./S1));

% calculate diameter of light beam at 1 m
dl = 2*1*tand(-theta2);

figure; plot(S1*1e3,theta2); xlabel('S1 [mm]'); ylabel('\theta_2 [deg]'); grid on
figure; plot(S1*1e3,dl); xlabel('S1 [mm]'); ylabel('d_l [m]'); grid on

% fprintf('S1 = %1.3f m, theta2 = %3.1f deg\n', S1, theta2)
