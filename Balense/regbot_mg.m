%% Simscape multibody model og Regbot in balance
% initial setup with motor velocity controller 
% this is intended as simulation base for balance control.
%
close all
clear

%%
function gamma_M = overshoot2phase_margin(M_p)
% Returns phase margin in degrees given overshoot (not in percentage)
% M_p must be positive
    zeta = log(1/M_p)/sqrt(pi^2 + log(1/M_p)^2);
    gamma_M = rad2deg(atan(2*zeta/sqrt(-2*zeta^2 + sqrt(1 + 4*zeta^4))));
end

function zeta = overshoot2damping(M_p)
% Returns damping ratio 0 < zeta < 1 for Overshoot M_p (not in percentage)
    zeta = log(1/M_p)/sqrt(pi^2 + log(1/M_p)^2);
end
function omega_c = crossover_second_order(omega_n,zeta)
% Returns the crossover frequency of a second order system given its
% natural frequency and its damping ratio
    omega_c = omega_n*sqrt(-2*zeta^2 + sqrt(1 + 4*zeta^4));
end
%% Simulink model name
model='regbot_1mg';

%% parameters for REGBOT
% motor
RA = 3.3/2;    % ohm (2 motors)
JA = 1.3e-6*2; % motor inertia
LA = 6.6e-3/2; % rotor inductor (2 motors)
BA = 3e-6*2;   % rotor friction
Kemf = 0.0105; % motor constant
Km = Kemf;
% køretøj
NG = 9.69; % gear
WR = 0.03; % wheel radius
Bw = 0.155; % wheel distance
% 
% model parts used in Simulink
mmotor = 0.193;   % total mass of motor and gear [kg]
mframe = 0.32;    % total mass of frame and base print [kg]
mtopextra = 0.97 - mframe - mmotor; % extra mass on top (charger and battery) [kg]
mpdist =  0.10;   % distance to lit [m]
% disturbance position (Z)
pushDist = 0.1; % relative to motor axle [m]

%% wheel velocity controller (no balance) PI-regulator
% sample (usable) controller values
Kpwv = 15;     % Kp
tiwv = 0.05;   % Tau_i
Kffwv = 0;     % feed forward constant
startAngle = 10;  % tilt in degrees at time zero
twvlp = 0.005;    % velocity noise low pass filter time constant (recommended)
Inputangel=30
Inputangelred=deg2rad(Inputangel)

%% Estimate transfer function for base system using LINEARIZE
% Motor volatge to wheel velocity (wv)
load_system(model);
open_system(model);
% define points in model
ios(1) = linio(strcat(model,'/vel_ref'),1,'openinput');
ios(2) = linio(strcat(model, '/tilt_angle'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model   
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gwv = minreal(tf(num, den))

%% Bodeplot
h = figure(100)
bode(Gwv)
grid on
title('Transfer function from motor voltage to velocity')
saveas(h, 'motor to velocity.png');
%%
t1 = 1/9.78
Ki=tf([t1,1],[t1,0])
Go=-Ki*Gwv
Gc=Go/(1+Go)

%%

%second order controler
%oversudte 20%
Mp=0.05
%rise time can også gøres med seteling time så skal wn bare laves om til
ts=0.5
Wn=4/(Cc*ts) %hvor ts er seteling time
Ni = 3

%%tr=2
%Wn=1.8/tr
%%

%
Cc=overshoot2damping(Mp)
Wc=crossover_second_order(Wn,Cc)
Ym=overshoot2phase_margin(Mp)
%%
G_eval = evalfr(Go, 1j*Wc);

theta = angle(G_eval); % Phase in radians in wc Rad/s
magnitude = abs(G_eval); % Get the magnitude in wc Rad/s

alpha_sol = -(sin(-Ym - atan(1/3) + theta) - 1)/(1 + sin(-Ym - atan(1/3) + theta))
%%
tau_d = 1/(sqrt(alpha_sol)*Wc)
C_lead = tf([tau_d, 1], [alpha_sol*tau_d, 1])


Kp=1/magnitude

tau_i=Ni/Wc

C_i=tf([tau_i,1],[tau_i,0])

Cpil = Kp*C_PI * C_lead
