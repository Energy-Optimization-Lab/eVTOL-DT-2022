%% Emrax_params.m
% By Jack Hannum, University of South Carolina 
% April 2023
% Based on Imperix TN111, reparameterized for the
% Emrax 228 MV CC Permanent Magnet Synchronous Machine
% (PMSM)
close all; clc;
fprintf('Emrax_params.m');
%% =================== System parameters ======================== %%
ControlFreq    = 20e3;           % Control interrupt frequency [Hz]
ControlPeriod  = 1/ControlFreq;  % Control interrupt period [s]
SamplingPeriod = ControlPeriod;  % Sampling period [s]

Ipmsm.Ts = ControlPeriod;        % Current controller control period [s]

% Deadtime
% dt = 1.2*(13e-9 + 20e-9 + 155e-9 + 380e-9);
dt = 1e-6;
%% =================== Plant model ======================== %%

% PMSM
Pmsm.p = 10;			% pole pairs [-]
Pmsm.Usn = 470;			% rated voltage [V]
Pmsm.Isn = 160;         % rated current [Arms]
Pmsm.fsn = 83.3;	    % rated frequency [Hz] (5000 RPM)
Pmsm.wsn = 2*pi*83.3;	% rated angular frequency [rad/s]
% Pmsm.Rs = 8e-3;       % stator resistance [Ohm]
% Pmsm.Rs = 0.06;       % stator resistance [Ohm]
Pmsm.Rs = 1.1;         % stator resistance [Ohm]
% Pmsm.Ld = 0.5*2.5e-3;        % d-axis inductance [H]
% Pmsm.Lq = 0.5*2.5e-3;        % q-axis inductance [H]
Pmsm.Lext = 2.5e-3;
Pmsm.Ld = Pmsm.Lext;        % d-axis inductance [H]
Pmsm.Lq = Pmsm.Lext;        % q-axis inductance [H]
Pmsm.Jm = 0.0421;		% moment of inertia [Nm^2]
Pmsm.kF = 0;            % friction coefficient 
Rl = 1.2-2*8e-3;
% Pmsm.Phi_pm = 0.75/sqrt(2);   % flux from permanent magnets [Nm/Apeak]
Pmsm.Phi_pm = 0.0355;

Pmsm.Temn = 13; 		% nominal torque [Nm]
Pmsm.Pn = 42e3;         % nominal active power [W]

Pmsm.Is_max = 340;		% [A] peak current
Pmsm.Wm_Max = 6500;		% [rpm]

% Encoder parameters
PPR = 8192/4;             % Encoder pulses per turn [-]

Vdc0 = 400;      % [V]
%% =============== PMSM - Current control ================= %%
s = tf('s');
pade = (1-s*ControlPeriod/4)/(1+s*ControlPeriod/4);
% Hd = (1/(Pmsm.Rs+Rl))/(1+s*Pmsm.Ld/(Pmsm.Rs+Rl));
% Hq = (1/(Pmsm.Rs+Rl))/(1+s*Pmsm.Lq/(Pmsm.Rs+Rl));
Hd = pade*(1/(Pmsm.Rs))/(1+s*Pmsm.Ld/(Pmsm.Rs));
Hq = pade*(1/(Pmsm.Rs))/(1+s*Pmsm.Lq/(Pmsm.Rs));
% Hd = c2d(Hd,ControlPeriod,'tustin');
% Hq = c2d(Hq,ControlPeriod,'tustin');
f_i = ControlFreq/10;
opts = pidtuneOptions('PhaseMargin',60);
PI_d = pidtune(Hd,'PI',2*pi*f_i,opts)
PI_q = pidtune(Hq,'PI',2*pi*f_i,opts)
% Plant model
Ipmsm.K1 = 1 / Pmsm.Rs;
Ipmsm.T1 = Pmsm.Ld / Pmsm.Rs;
Ipmsm.K2 = 1 / Pmsm.Rs;
Ipmsm.T2 = Pmsm.Lq / Pmsm.Rs;

% Small time constants
Ipmsm.Tctrl = Ipmsm.Ts;			% Control delay [s]
Ipmsm.Tpwm = 0.5 * ControlPeriod;		% Modulator delay [s]
Ipmsm.Tsens = 0;                        % No filter used, ideal sensor
% Ipmsm.Tsens = 1 / (2*pi*1.6e3);		% Sensing delay (if LPF on ADC) [s]

Ipmsm.Kcm = 1;
% Ipmsm.Ttot = Ipmsm.Tctrl + Ipmsm.Tpwm + Ipmsm.Tsens;
Ipmsm.Ttot = ControlPeriod + 0.5*ControlPeriod;
% Magnitude optimum
Ipmsm.Tn_d = Ipmsm.T1;
Ipmsm.Ti_d = 2 * Ipmsm.Kcm * Ipmsm.K1 * Ipmsm.Ttot;
Ipmsm.Tn_q = Ipmsm.T2;
Ipmsm.Ti_q = 2 * Ipmsm.Kcm * Ipmsm.K2 * Ipmsm.Ttot;

Ipmsm.Kp_d = Ipmsm.Tn_d / Ipmsm.Ti_d;
Ipmsm.Ki_d = 1 / Ipmsm.Ti_d;
Ipmsm.Kp_q = Ipmsm.Tn_q / Ipmsm.Ti_q;
Ipmsm.Ki_q = 1 / Ipmsm.Ti_q;

%% Speed Controller
s = tf('s');
fc_s = 5;
Hs = 1/(s*Pmsm.Jm);
% Hs = c2d(Hs,ControlPeriod,'tustin');
opts = pidtuneOptions('PhaseMargin',80);
PI_s = pidtune(Hs,'PI',2*pi*fc_s,opts)
T = PI_s*Hs;
figure;margin(T);
HCL = feedback(T,1);
figure;step(188*HCL);

%% SVPWM Initialization
Cclarke = sqrt(2/3)*[1,-0.5,-0.5;0, sqrt(3)/2, sqrt(3)/2; 1/sqrt(2),1/sqrt(2),1/sqrt(2);];
Iinit = 5;
Rs = 2*Pmsm.Rs;
Va = Iinit*Rs;
Vabc = [Va;0;0];
Vab0 = Cclarke*Vabc;
Valpha = Vab0(1)
Vbeta = Vab0(2)
Vzero = Vab0(3)
Kb = 0.1;
