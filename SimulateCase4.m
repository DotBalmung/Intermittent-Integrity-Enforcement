clear all;close all;clc;

%% define parameters for the simulation
rng(777)
tau = 2.5;%s - time gap distance
Tl = 0.45; %eq (1), related to generalized vehicle longitudinal dynamic
Kl = 1;
Ac = [0 1 -tau ; 0 0 -1 ; 0 0 -1/Tl];
Bc = [0 0 ; 0 1 ; Kl/Tl 0];
Cc = eye(3);
Dc = 0;

%%  Discretize the model with sampling rate Ts
Ts = 0.02;
csys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(csys,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%% rest of parameters
Gatk = eye(3); % Matrix that maps attacks to corresponding positions
W = eye(3); % Used for LQR controler cost function
U = eye(2); % Used for LQR controler cost fucntion
Q = 0.01*eye(3); % Covariance matrix of state noise
R = 0.01*eye(3); % Covariance matrix of sensor noise
AtkDuration=1000;
x0 = [0;0;0];
IsAttk = 1; % Specifies whether or not to apply attack
% Since the attack searching algorithm is finding suboptimal attacks, we
% need to scale obtained attack vectors to emulate optimal attacker (i.e.,
% the attacker that reaches e_max in every step of the attack). For
% different values of L, we will have different scaling factor.
%AttkScale=0.45/1.3911; %i'm supposed to hit 0.45 error in position for L=2
%AttkScale=0.56/1.8655; %i'm supposed to hit 0.56 error in position for L=3
%AttkScale=0.61/2.0366; %i'm supposed to hit 0.61 error in position for L=4
AttkScale=0.67/2.1972; %i'm supposed to hit 0.67 error in position for L=5
%AttkScale=0.65/2.4238; %i'm supposed to hit 0.72 error in position for L=6
%AttkScale=1/4.5421; %i'm supposed to hit 1 error in position for L=13
refx = zeros(3,2000); % reference value for the controller
authT = 5; %How often the data is being authenticated, same as L from above

[x_ea,y_ea,ya_ea,x_hat_ea,u_ea,a_ea] = SimulateClosedLoopSystemV2GAasEA(A,B,C,Gatk,W,U,Q,R,x0,IsAttk,AtkDuration,AttkScale,refx,authT);
beep
