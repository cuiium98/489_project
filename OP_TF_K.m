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
K = 5.12;

s = tf('s');

% =========================
% Transfer functionsi 
% =========================
Gc = 1;

Glc = 1 / (L*C*s^2 + (L/R)*s + 1);

% CORRECTED: PWM modulator gain (continuous-time, no 1/2fsw scaling)
% Half-period transport delay approximated as Pade 1st order
Td     = 1 / (2*fsw);                          % half-period delay = 2.5us
Gpwm   = (Vc / (2*Vtri)) * (1 - s*Td/2) / (1 + s*Td/2);  % Pade approx of e^(-sTd)

H = Rfb2 / (Rfb1 + Rfb2);

T = K*(Gc * Glc * Gpwm * H);

% =========================
% Bode plot
% =========================
figure;
bode(T);
grid on;
title('Bode Plot of T(s) — Added Gain K');

% =========================
% Nyquist plot
% =========================
figure;
nyquist(T);
grid on;
title('Nyquist Plot of T(s) — Added Gain K');

% =========================
% Root locus
% =========================
figure;
rlocus(T);
grid on;
title('Root Locus of T(s) — Added Gain K');

% =========================
% Print gain and phase margins
% =========================
[Gm, Pm, Wcg, Wcp] = margin(T);
fprintf('Gain Margin:       %.2f dB  at %.2f krad/s\n', 20*log10(Gm), Wcg/1e3);
fprintf('Phase Margin:      %.2f deg  at %.2f krad/s\n', Pm, Wcp/1e3);
fprintf('DC Loop Gain T(0): %.4f (%.2f dB)\n', abs(dcgain(T)), 20*log10(abs(dcgain(T))));