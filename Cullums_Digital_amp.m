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
Amp = (Vbus/2*Vtripk);      % Modulator gain

s = tf('s');
G_lc = 1/(L1*C1*s^2 + (L1/Rload)*s + 1);        % 2nd-order LC lowpass with load damping

% Plant TF

G_p = Amp * G_lc * Afb;

%% Discretize Plant (ZOH)
Ts = 1/fsw;                      % sample period, s
G_pd = c2d(G_p, Ts, 'zoh');   % discretize plant, ZOH captures sample-hold exactly


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
Gc_shape = (1 + s/wz)^2 / (s*(1 + s/wp)^2);         % compensator shape (gain set below)
magc = squeeze(bode(Gc_shape*G_p, wc));           % loop magnitude at wc before scaling
Kgain = 1/magc;                                   % gain that forces |loop| = 1 at wc
Gc = Kgain * Gc_shape;                            % final Type III compensator

%% Discretize Compensator (Tustin, prewarped at crossover)
opts = c2dOptions('Method', 'tustin', 'PrewarpFrequency', wc);
G_cd = c2d(Gc, Ts, opts); % G_comp = your continuous Type III / PR TF


%% Discrete PR Controller (Tustin) — coefficients per figure
f2k = 2000;
Tsd = 1/(2*fsw);        % sample period, s
Ki = 170.4;        % resonant (integral) gain
Kp = 0;        % proportional gain
w0 = (2*pi*f2k);        % resonant frequency, rad/s   (2*pi*f0)
wc = (2*pi*10);        % resonant bandwidth,  rad/s

% Coefficients (exactly as in the figure)
a1 = 4*Ki*Ts*wc;
b0 = Ts^2*w0^2 + 4*Ts*wc + 4;
b1 = 2*Ts^2*w0^2 - 8;
b2 = Ts^2*w0^2 - 4*Ts*wc + 4;

%% Transfer-function form (for Bode / margin analysis)
z    = tf('z', Ts);
Gcr  = a1*(1 - z^-2) / (b0 + b1*z^-1 + b2*z^-2);   % resonant term
G_prd  = Kp + Gcr;  

%% ---------- Bode plot: Resonator
figure;
bode(G_prd);
legend('Resonator Bode Plot');
grid on;
title('Open-Loop: effect of resonant term at 2 kHz');% full PR controller

%% ---------- Closed-loop analysis (Type III only) ----------
L                = G_cd * G_pd;          % compensated open-loop
L_new                = G_cd * G_pd * G_prd;          % compensated open-loop
[GM, PM, ~, wcp] = margin(L);        % margins
GM_dB            = 20*log10(GM);     % gain margin             [dB]

fprintf('PM = %.1f deg\n', PM);
fprintf('GM = %.1f dB\n',  GM_dB);
%% ---------- Bode plot: Type III vs Type III + resonant ----------
figure;
bode(L, L_new);
legend('Type III only', 'Type III + resonant');
grid on;
title('Open-Loop: effect of resonant term at 2 kHz');











