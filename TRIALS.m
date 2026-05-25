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

% Plant/loop without compensator
Gp = Glc * Gpwm * H;

% =========================
% Desired crossover and phase margin
% =========================
fc = 20e3;              % [Hz]
wc = 2*pi*fc;           % [rad/s]

PM_desired = 60;        % [deg]

% =========================
% Uncompensated plant at wc
% =========================
resp_Gp = evalfr(Gp, 1j*wc);

Gp_mag = abs(resp_Gp);
Gp_gain_dB = 20*log10(Gp_mag);
phi_sys = angle(resp_Gp) * 180/pi;

% unwrap phase safely into negative range
while phi_sys > 0
    phi_sys = phi_sys - 360;
end

while phi_sys <= -360
    phi_sys = phi_sys + 360;
end

% =========================
% Required boost and n-factor
% =========================
% Type-III has an integrator, so it contributes -90 degrees.
% Required boost:
phi_boost_required = PM_desired - phi_sys - 90;

% Do not allow negative boost
if phi_boost_required < 0
    phi_boost_required = 0;
end

% Type-III max theoretical boost tends to 180 deg
if phi_boost_required >= 180
    error('Required phase boost is too large. Lower fc or change plant/compensator design.');
end

% For Type-III with two identical lead sections:
% phi_boost = 2*asin((n-1)/(n+1))
% where n = wp/wz

n_min = (1 + sind(phi_boost_required/2)) / ...
        (1 - sind(phi_boost_required/2));

% Clamp n_min to at least 1
n_min = max(n_min, 1);

% =========================
% Manual n selection
% =========================
% n is the spreading ratio:
% n = wp / wz
%
% Increase n to separate zeros and poles more.
% n = n_min gives the minimum required boost.

n = 1.25;

% Clamp n to at least 1
n = max(n, 1);

if n < n_min
    warning('Chosen n is less than n_min. Desired phase margin may not be achieved.');
end

% =========================
% Zero and pole locations
% =========================
% wc must be geometric mean:
% wc = sqrt(wz*wp)

wz = wc / sqrt(n);
wp = wc * sqrt(n);

fz = wz / (2*pi);
fp = wp / (2*pi);

% Verify geometric mean
assert(abs(sqrt(wz*wp) - wc)/wc < 1e-10, ...
       'Geometric mean check failed: wc is not sqrt(wz*wp).');

% Nyquist check
f_nyq = fsw/2;
w_nyq = 2*pi*f_nyq;

if wp >= w_nyq
    error('wp is above Nyquist. Reduce n or reduce fc. wp must be below fsw/2.');
end

% Actual phase boost from chosen n
phi_boost_actual = 2 * asind((n - 1) / (n + 1));

% =========================
% Type-III compensator with Kc = 1
% =========================
Gc_shape = ((1 + s/wz)^2) / (s * (1 + s/wp)^2);

Loop_Kc1 = Gc_shape * Gp;

resp_Loop_Kc1 = evalfr(Loop_Kc1, 1j*wc);

Loop_Kc1_mag = abs(resp_Loop_Kc1);
Loop_Kc1_gain_dB = 20*log10(Loop_Kc1_mag);
Loop_Kc1_phase = angle(resp_Loop_Kc1) * 180/pi;

while Loop_Kc1_phase > 0
    Loop_Kc1_phase = Loop_Kc1_phase - 360;
end

while Loop_Kc1_phase <= -360
    Loop_Kc1_phase = Loop_Kc1_phase + 360;
end

% Required gain to force crossover at wc
Kc = 1 / Loop_Kc1_mag;

% Final compensator and loop
Gc = Kc * Gc_shape;
T = Gc * Gp;

% =========================
% Final response at designed wc
% =========================
resp_T = evalfr(T, 1j*wc);

T_mag = abs(resp_T);
T_gain_dB = 20*log10(T_mag);
T_phase = angle(resp_T) * 180/pi;

while T_phase > 0
    T_phase = T_phase - 360;
end

while T_phase <= -360
    T_phase = T_phase + 360;
end

PM_expected = 180 + T_phase;

% =========================
% Slope estimate at wc
% =========================
% Integrator: -20 dB/dec
% Two lead sections contribution:
% +40*(n-1)/(n+1)

slope_Gc_wc = -20 + 40*((n - 1)/(n + 1));

% =========================
% Print values
% =========================
fprintf('\n========== Type-III Compensator with Manual n ==========\n');

fprintf('\nDesired Design Targets:\n');
fprintf('fc:                         %.2f kHz\n', fc/1e3);
fprintf('wc:                         %.4e rad/s\n', wc);
fprintf('Desired PM:                  %.2f deg\n', PM_desired);

fprintf('\nUncompensated Plant at wc:\n');
fprintf('|Gp(jwc)|:                   %.6f\n', Gp_mag);
fprintf('Gain:                        %.2f dB\n', Gp_gain_dB);
fprintf('Phase phi_sys:               %.2f deg\n', phi_sys);

fprintf('\nType-III n Design:\n');
fprintf('Required phi_boost:          %.2f deg\n', phi_boost_required);
fprintf('Minimum n from formula:      %.6f\n', n_min);
fprintf('Chosen manual n:             %.6f\n', n);
fprintf('Actual phi_boost from n:     %.2f deg\n', phi_boost_actual);
fprintf('Actual wp/wz ratio:          %.6f\n', wp/wz);

fprintf('\nZero/Pole Locations:\n');
fprintf('wz1 = wz2:                   %.4e rad/s\n', wz);
fprintf('wp1 = wp2:                   %.4e rad/s\n', wp);
fprintf('fz1 = fz2:                   %.2f Hz\n', fz);
fprintf('fp1 = fp2:                   %.2f Hz\n', fp);
fprintf('Geometric mean frequency:    %.2f Hz\n', sqrt(fz*fp));
fprintf('Nyquist frequency:           %.2f Hz\n', f_nyq);

fprintf('\nKc Calculation:\n');
fprintf('|Loop with Kc=1| at wc:      %.6e\n', Loop_Kc1_mag);
fprintf('Gain with Kc=1 at wc:        %.2f dB\n', Loop_Kc1_gain_dB);
fprintf('Phase with Kc=1:             %.2f deg\n', Loop_Kc1_phase);
fprintf('Required Kc = 1/mag:         %.6e\n', Kc);

fprintf('\nFinal Compensated Loop at Designed wc:\n');
fprintf('|T(jwc)|:                    %.6f\n', T_mag);
fprintf('Gain:                        %.2f dB\n', T_gain_dB);
fprintf('Phase:                       %.2f deg\n', T_phase);
fprintf('Expected PM at designed fc:  %.2f deg\n', PM_expected);

fprintf('\nCompensator Shape:\n');
fprintf('Approx Gc slope at wc:       %.2f dB/dec\n', slope_Gc_wc);

if slope_Gc_wc > 0
    fprintf('Gc magnitude is rising at wc.\n');
elseif slope_Gc_wc < 0
    fprintf('Gc magnitude is falling at wc.\n');
else
    fprintf('Gc magnitude is flat at wc.\n');
end

% =========================
% Margins
% =========================
[Gm, Pm, Wcg, Wcp] = margin(T);

fprintf('\nMATLAB Critical Margin Results:\n');

if isinf(Gm)
    fprintf('Gain Margin:                 Inf\n');
else
    fprintf('Gain Margin:                 %.2f dB at %.2f krad/s\n', 20*log10(Gm), Wcg/1e3);
end

fprintf('Critical Phase Margin:       %.2f deg at %.2f krad/s\n', Pm, Wcp/1e3);
fprintf('Critical Gain Crossover:     %.2f kHz\n', Wcp/(2*pi*1e3));

fprintf('\nDesigned Crossover Check:\n');
fprintf('Designed Gain Crossover:     %.2f kHz\n', fc/1e3);
fprintf('|T(jwc)| at designed fc:     %.6f\n', T_mag);
fprintf('Gain at designed fc:         %.2f dB\n', T_gain_dB);

% =========================
% Closed-loop response and ESS
% =========================
T_cl = feedback(T, 1);

info = stepinfo(T_cl);
final_value = dcgain(T_cl);
ess = abs(1 - final_value);

fprintf('\nClosed-Loop Step Response:\n');
fprintf('Final Value:                 %.6f\n', final_value);
fprintf('ESS:                         %.6e\n', ess);
fprintf('Rise Time:                   %.6e s\n', info.RiseTime);
fprintf('Settling Time:               %.6e s\n', info.SettlingTime);
fprintf('Overshoot:                   %.2f %%\n', info.Overshoot);
fprintf('Peak Value:                  %.6f\n', info.Peak);
fprintf('Peak Time:                   %.6e s\n', info.PeakTime);

% =========================
% Plots
% =========================
figure;
bode(Gp);
grid on;
title('Uncompensated Plant Gp(s)');

figure;
bode(Gc);
grid on;
title('Type-III Compensator Only Gc(s)');

figure;
bode(T);
grid on;
title('Compensated Open-Loop T(s)');

figure;
margin(T);
grid on;
title('Margins of Compensated Loop');

figure;
step(T_cl);
grid on;
title('Closed-Loop Step Response');

figure;
nyquist(T);
grid on;
title('Nyquist Plot');

figure;
rlocus(T);
grid on;
title('Root Locus');