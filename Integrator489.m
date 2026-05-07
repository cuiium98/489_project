clc;
clear;
close all;

% =========================
% Parameters
% =========================
Vc   = 5;          % [V]
Vtri = 1.6;        % [V]
fsw  = 200e3;      % [Hz]
L    = 22.5e-6;    % [H]
C    = 0.7e-6;     % [F]
R    = 8;          % [ohm]
Rfb1 = 10e3;       % [ohm]
Rfb2 = 1.429e3;    % [ohm]

s = tf('s');

% =========================
% Transfer functions
% =========================
Gc   = 1;

Glc  = 1 / (L*C*s^2 + (L/R)*s + 1);

Gpwm = (Vc / (2*Vtri)) * (1 / (2*fsw));

H    = Rfb2 / (Rfb1 + Rfb2);

T    = Gc * Glc * Gpwm * H;



% =========================
% Bode plot
% =========================
figure;
bode(T);
grid on;
title('Bode Plot of T(s)');

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