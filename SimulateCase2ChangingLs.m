clear all;close all;clc;

%% define parameters for the simulation
rng(777)
m = 1573; %kg
Iz = 2873; %vehicle yaw moment of inertia
lf = 1.1; %length in front of c.g.
lr = 1.58; % length behind of c.g. (rear)
Caf = 80000; %Cornering coefficient front
Car = 80000; %Cornering coefficient rear
Vx = 20; %mps
Ac = [0 1 0 0;0 -2*(Caf+Car)/(m*Vx) 2*(Caf+Car)/m 2*(-Caf*lf+Car*lr)/(m*Vx); 0 0 0 1 ; 0 -2*(Caf*lf-Car*lr)/(Iz*Vx) 2*(Caf*lf-Car*lr)/Iz -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
Bc = [0;2*Caf/m;0;2*Caf*lf/Iz];
Cc = eye(4);
Dc = 0;

addProb = 0.0001;
FalsePosProb = 0.05;

%%  Discretize the model with sampling rate Ts
Ts = 0.02;
csys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(csys,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%% rest of parameters
Gatk = eye(4); % Matrix that maps attacks to corresponding positions
W = eye(4); % Used for LQR controler cost function
U = 1; % Used for LQR controler cost fucntion
Q = 0.001*eye(4); % Covariance matrix of state noise
R = 0.001*eye(4); % Covariance matrix of sensor noise
AtkDuration=1000;
x0 = [0;0;0;0];
IsAttk = 1; % Specifies whether or not to apply attack
AttkScale = 1; % Kept at 1, as we are using attack patterns here
refx = zeros(4,2000); % reference value for the controller

load('LPlistOfLs.mat'); % Loads the list of L values computed by the scheduling scripts
[x_ea,y_ea,ya_ea,x_hat_ea,u_ea,a_ea] = SimulateClosedLoopSystemChangingL(A,B,C,Gatk,W,U,Q,R,x0,IsAttk,AtkDuration,AttkScale,refx,LanePositionLsList,0);
beep
