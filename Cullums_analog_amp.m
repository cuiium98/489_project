%% Hardware Platform Parameters

clear; clc;

% Modulator Stage
Vtripk = 0.9;

% Output filter - Stage 1 (in-loop)
L1 = 33e-6;        % Stage 1 inductance, H
C1 = 500e-9;        % Stage 1 capacitance, F

% Power stage
Vbus = 12;      % DC bus voltage, V
fsw  = 200e+3;      % Switching frequency, Hz
Rload = 8;     % Load resistance, Ohm

% Feedback
Afb = 0.233;         % Feedback divider gain

% Gate drive / control
Td = (1/(2*fsw));              % Single-Rate Update delay
[num,den] = pade(Td, 2);     % 2nd-order Pade approximation of e^(-sTd)
G_d = tf(num,den);           % rational delay model
Amp = (Vbus/2*Vtripk);      % Modulator gain

s = tf('s');
G_lc = 1/(L1*C1*s^2 + (L1/Rload)*s + 1);        % 2nd-order LC lowpass with load damping

% Plant TF

G_p = Amp * G_d * G_lc * Afb;

%% Type III Compensator - K-Factor Method

boost = 23;     % required phase boost at fc, degrees
fc    = 20e+3;     % desired crossover frequency, Hz
wc = 2*pi*fc;


K  = (tand(boost/4 + 45))^2;   % K-factor
wz = wc / sqrt(K);              % coincident zero frequency, Hz
wp = wc * sqrt(K);               % coincident pole frequency, Hz

fprintf('Compensator Pole location: %f rads\n', wp);
fprintf('Compensator Zero location: %f rads\n', wz);

%% --- Build the Type III compensator (integrator + zero + pole) ---
Gc_shape = (1 + s/wz)^2 / (s*(1 + s/wp)^2);          % compensator shape (gain set below)
magc = squeeze(bode(Gc_shape*G_p, wc));           % loop magnitude at wc before scaling
Kgain = 1/magc;                                   % gain that forces |loop| = 1 at wc
Gc = Kgain * Gc_shape;                            % final Type III compensator

%% ---------- Closed-loop analysis (Type III only) ----------
L                = Gc * G_p;          % compensated open-loop
[GM, PM, ~, wcp] = margin(L);        % margins
GM_dB            = 20*log10(GM);     % gain margin             [dB]

fprintf('PM = %.1f deg\n', PM);
fprintf('GM = %.1f dB\n',  GM_dB);
%% Phase of System at a Specific Frequency

sys = G_p;       % transfer function object
f   = 20e+3;       % frequency of interest, Hz

w = 2*pi*f;
[~, phase] = bode(sys, w);
phase = squeeze(phase);

fprintf('Pre compensator Phase at 20kHz: %f deg\n', phase);











