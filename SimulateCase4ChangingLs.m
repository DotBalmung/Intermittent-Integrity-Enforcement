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
AttkScale=1; % Kept at 1, as we are using attack patterns here
refx = zeros(3,2100); % reference value for the controller

load('CClistOfLs.mat'); % Loads the list of L values computed by the scheduling scripts
[x_ea,y_ea,ya_ea,x_hat_ea,u_ea,a_ea] = SimulateClosedLoopSystemChangingL(A,B,C,Gatk,W,U,Q,R,x0,IsAttk,AtkDuration,AttkScale,refx,CruiseControlLsList,1);
beep
