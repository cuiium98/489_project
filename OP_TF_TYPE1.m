clc; clear; close all;

% =========================
% Parameters
% =========================
Vc   = 5;        
Vtri = 1.6;      
fsw  = 200e3;    
L    = 22.5e-6;  
C    = 0.7e-6;   
R    = 8;        
Rfb1 = 10e3;     
Rfb2 = 1.429e3;  

s = tf('s');

% =========================
% Plant transfer functions
% =========================
Glc = 1 / (L*C*s^2 + (L/R)*s + 1);

Td = 1 / (2*fsw);
Gpwm = (Vc / (2*Vtri)) * (1 - s*Td/2) / (1 + s*Td/2);

H = Rfb2 / (Rfb1 + Rfb2);

% Plant seen by compensator
Gp = Glc * Gpwm * H;

% =========================
% Desired crossover
% =========================
fc = 20e3;              % [Hz]
wc = 2*pi*fc;           % [rad/s]

% Compute |Gp(jwc)|
Gp_mag = abs(freqresp(Gp, wc));

% For Gc = K/s:
% |(K/jwc)Gp(jwc)| = 1
K = wc / Gp_mag;

% Compensator
Gc = K / s;

% Final open-loop transfer function
T = Gc * Gp;

% =========================
% Closed-loop transient response + ESS
% =========================
T_cl = feedback(Gc*Gp, 1);   % unity negative feedback closed-loop

figure;
step(T_cl);
grid on;
title('Closed-Loop Step Response');

info = stepinfo(T_cl);

fprintf('\nTransient Response Info:\n');
fprintf('Rise Time:        %.6e s\n', info.RiseTime);
fprintf('Settling Time:    %.6e s\n', info.SettlingTime);
fprintf('Overshoot:        %.2f %%\n', info.Overshoot);
fprintf('Peak Time:        %.6e s\n', info.PeakTime);
fprintf('Peak Value:       %.6f\n', info.Peak);

% Steady-state error for unit step input
ess = abs(1 - dcgain(T_cl));

fprintf('\nSteady-State Error:\n');
fprintf('ESS for unit step input: %.6e\n', ess);

% =========================
% Print results
% =========================
fprintf('Desired crossover frequency: %.2f kHz\n', fc/1e3);
fprintf('Desired crossover frequency: %.2f krad/s\n', wc/1e3);
fprintf('|Gp(jwc)|: %.6f\n', Gp_mag);
fprintf('Required K: %.6e\n', K);

[Gm, Pm, Wcg, Wcp] = margin(T);

fprintf('\nAfter applying Gc = K/s:\n');
fprintf('Gain Margin:       %.2f dB at %.2f krad/s\n', 20*log10(Gm), Wcg/1e3);
fprintf('Phase Margin:      %.2f deg at %.2f krad/s\n', Pm, Wcp/1e3);
fprintf('Actual crossover:  %.2f kHz\n', Wcp/(2*pi*1e3));

% =========================
% Plots
% =========================
figure;
bode(T);
grid on;
title('Open-Loop Bode Plot with Gc = K/s');

figure;
nyquist(T);
grid on;
title('Nyquist Plot with Gc = K/s');

figure;
rlocus(T);
grid on;
title('Root Locus with Gc = K/s');