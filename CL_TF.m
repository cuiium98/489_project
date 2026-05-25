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
Rfb1 = 10e3;     % [ohm]
Rfb2 = 1.429e3;  % [ohm]
s = tf('s');
% =========================
% Transfer functions
% =========================
Gc = 1;
Glc = 1 / (L*C*s^2 + (L/R)*s + 1);
Td     = 1 / (2*fsw);
Gpwm   = (Vc / (2*Vtri)) * (1 - s*Td/2) / (1 + s*Td/2);
H = 5*(Rfb2 / (Rfb1 + Rfb2));
T = Gc * Glc * Gpwm * H;

% =========================
% Closed Loop
% =========================
CL = feedback(T, 1);

% =========================
% Steady State Error
% =========================
dc_gain = dcgain(CL);
sse = 1 - dc_gain;
fprintf('Closed Loop DC Gain:  %.6f\n', dc_gain);
fprintf('Steady State Error:   %.6f (%.4f%%)\n', sse, sse*100);

% =========================
% Step Response
% =========================
figure;
step(CL);
grid on;
title('Closed Loop Step Response');

% =========================
% Bode plot
% =========================
figure;
bode(T);
grid on;
title('Bode Plot of T(s) — Open Loop');

% =========================
% Nyquist plot
% =========================
figure;
nyquist(T);
grid on;
title('Nyquist Plot of T(s)');

% =========================
% Root locus
% =========================
figure;
rlocus(T);
grid on;
title('Root Locus of T(s)');

% =========================
% Print gain and phase margins
% =========================
[Gm, Pm, Wcg, Wcp] = margin(T);
fprintf('Gain Margin:       %.2f dB  at %.2f krad/s\n', 20*log10(Gm), Wcg/1e3);
fprintf('Phase Margin:      %.2f deg  at %.2f krad/s\n', Pm, Wcp/1e3);
fprintf('DC Loop Gain T(0): %.4f (%.2f dB)\n', abs(dcgain(T)), 20*log10(abs(dcgain(T))));