clc; clear; close all;

% =========================
% Parameters
% =========================
Vc   = 5;        % [V]
Vtri = 1.6;      % [V]
fsw  = 200e3;    % [Hz]
L    = 22.5e-6;  % [H]
C    = 0.7e-6;   % [F]
R    = 8;        % [ohm]
Rfb1 = 10e3;     
Rfb2 = 1.429e3;  

s = tf('s');

% =========================
% Plant transfer functions
% =========================
Glc = 1 / (L*C*s^2 + (L/R)*s + 1);

Td = 1 / (2*fsw);

Gpwm = (Vc / (2*Vtri)) * ...
       (1 - s*Td/2) / (1 + s*Td/2);

H = Rfb2 / (Rfb1 + Rfb2);

% Uncompensated plant/loop
Gp = Glc * Gpwm * H;

% =========================
% Desired crossover and phase margin
% =========================
fc = 20e3;              % [Hz]
wc = 2*pi*fc;           % [rad/s]

PM_desired = 60;        % [deg] choose 45 to 60 usually

% =========================
% Plant gain and phase at wc
% =========================
resp_Gp = freqresp(Gp, wc);

Gp_mag = abs(resp_Gp);
Gp_gain_dB = 20*log10(Gp_mag);
phi_sys = angle(resp_Gp) * 180/pi;

% =========================
% Type-II compensator design
% =========================
phi_boost = PM_desired - phi_sys - 90;

k = tan(deg2rad(phi_boost/2 + 45));

wz = wc / k;
wp = k * wc;

fz = wz / (2*pi);
fp = wp / (2*pi);

% Type-II compensator without gain Kc
Gc_shape = (1 + s/wz) / (s * (1 + s/wp));

% Loop gain with Kc initially set to 1
Loop_Kc1 = Gc_shape * Gp;

resp_Loop_Kc1 = freqresp(Loop_Kc1, wc);
Loop_Kc1_mag = abs(resp_Loop_Kc1);
Loop_Kc1_gain_dB = 20*log10(Loop_Kc1_mag);
Loop_Kc1_phase = angle(resp_Loop_Kc1) * 180/pi;

% Required gain Kc to force crossover at wc
Kc = 1 / Loop_Kc1_mag;

% Final Type-II compensator
Gc = Kc * Gc_shape;

% Final compensated open-loop transfer function
T = Gc * Gp;

% =========================
% Final loop response at wc
% =========================
resp_T = freqresp(T, wc);

T_mag = abs(resp_T);
T_gain_dB = 20*log10(T_mag);
T_phase = angle(resp_T) * 180/pi;

PM_expected = 180 + T_phase;

% =========================
% Print all calculated values
% =========================
fprintf('\n========== Type-II Compensator Design ==========\n');

fprintf('\nDesired Design Targets:\n');
fprintf('fc:                 %.2f kHz\n', fc/1e3);
fprintf('wc:                 %.4e rad/s\n', wc);
fprintf('Desired PM:          %.2f deg\n', PM_desired);

fprintf('\nUncompensated Plant at wc:\n');
fprintf('|Gp(jwc)|:           %.6f\n', Gp_mag);
fprintf('Gain:                %.2f dB\n', Gp_gain_dB);
fprintf('Phase phi_sys:       %.2f deg\n', phi_sys);

fprintf('\nType-II Phase Boost Design:\n');
fprintf('Required phi_boost:  %.2f deg\n', phi_boost);
fprintf('k:                   %.6f\n', k);
fprintf('wz:                  %.4e rad/s\n', wz);
fprintf('wp:                  %.4e rad/s\n', wp);
fprintf('fz:                  %.2f Hz\n', fz);
fprintf('fp:                  %.2f Hz\n', fp);

fprintf('\nGain Setting:\n');
fprintf('|Loop with Kc=1|:    %.6e\n', Loop_Kc1_mag);
fprintf('Gain with Kc=1:      %.2f dB\n', Loop_Kc1_gain_dB);
fprintf('Phase with Kc=1:     %.2f deg\n', Loop_Kc1_phase);
fprintf('Required Kc:         %.6e\n', Kc);

fprintf('\nFinal Compensated Loop at wc:\n');
fprintf('|T(jwc)|:            %.6f\n', T_mag);
fprintf('Gain:                %.2f dB\n', T_gain_dB);
fprintf('Phase:               %.2f deg\n', T_phase);
fprintf('Expected PM:         %.2f deg\n', PM_expected);

% =========================
% Margins
% =========================
[Gm, Pm, Wcg, Wcp] = margin(T);

fprintf('\nMATLAB Margin Results:\n');
fprintf('Gain Margin:         %.2f dB at %.2f krad/s\n', 20*log10(Gm), Wcg/1e3);
fprintf('Phase Margin:        %.2f deg at %.2f krad/s\n', Pm, Wcp/1e3);
fprintf('Gain Crossover:      %.2f kHz\n', Wcp/(2*pi*1e3));

% =========================
% Closed-loop response and ESS
% =========================
T_cl = feedback(T, 1);

info = stepinfo(T_cl);
final_value = dcgain(T_cl);
ess = abs(1 - final_value);

fprintf('\nClosed-Loop Step Response:\n');
fprintf('Final Value:         %.6f\n', final_value);
fprintf('ESS:                 %.6e\n', ess);
fprintf('Rise Time:           %.6e s\n', info.RiseTime);
fprintf('Settling Time:       %.6e s\n', info.SettlingTime);
fprintf('Overshoot:           %.2f %%\n', info.Overshoot);
fprintf('Peak Value:          %.6f\n', info.Peak);
fprintf('Peak Time:           %.6e s\n', info.PeakTime);

% =========================
% Plots
% =========================
figure;
bode(Gp);
grid on;
title('Uncompensated Plant Gp(s)');

figure;
bode(T);
grid on;
title('Compensated Open-Loop T(s) with Type-II Compensator');

figure;
margin(T);
grid on;
title('Margins of Compensated Loop');

figure;
step(T_cl);
grid on;
title('Closed-Loop Step Response with Type-II Compensator');

figure;
nyquist(T);
grid on;
title('Nyquist Plot with Type-II Compensator');

figure;
rlocus(T);
grid on;
title('Root Locus with Type-II Compensator');

figure;
bode(Gc);
grid on;
title('Type-II Compensator Only');