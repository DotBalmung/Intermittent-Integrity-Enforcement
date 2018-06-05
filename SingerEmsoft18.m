clear all;close all;clc;

%% define parameters for the simulation
rng(777)
% m = 1573; %kg
% Iz = 2873; %vehicle yaw moment of inertia
% lf = 1.1; %length in front of c.g.
% lr = 1.58; % length behind of c.g. (rear)
% Caf = 80000; %Cornering coefficient front
% Car = 80000; %Cornering coefficient rear
% Vx = 20; %mps
% Ac = [0 1 0 0;0 -2*(Caf+Car)/(m*Vx) 2*(Caf+Car)/m 2*(-Caf*lf+Car*lr)/(m*Vx); 0 0 0 1 ; 0 -2*(Caf*lf-Car*lr)/(Iz*Vx) 2*(Caf*lf-Car*lr)/Iz -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
% Bc = [0;2*Caf/m;0;2*Caf*lf/Iz];
% Cc = eye(4);
% Dc = 0;

tau = 0.8;
Ac = [0 -1 0 ; 0 0 1 ; 0 0 -1/tau];
Bc = [0 ; 0 ; 1];
Cc = [1 0 0 ; 0 1 0];
Dc = 0;

addProb = 0.01;
FalsePosProb = 0.05;

%% Discretize the model with sampling rate Ts
Ts = 0.02;
csys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(csys,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%% rest of parameters
Gatk = eye(2); % Matrix that maps attacks to corresponding positions
W = eye(3); % Used for LQR controler cost function
U = 1; % Used for LQR controler cost fucntion
Q = 0.001*eye(3); % Covariance matrix of state noise
R = 0.001*eye(2); % Covariance matrix of sensor noise
AtkDuration=100;
x0 = [0;0;0];
IsAttk = 1; % Specifies whether or not to apply attack
% Since the attack searching algorithm is finding suboptimal attacks, we
% need to scale obtained attack vectors to emulate optimal attacker (i.e.,
% the attacker that reaches e_max in every step of the attack). For
% different values of L, we will have different scaling factor.
%AttkScale=0.047/0.8736; %i'm supposed to hit 0.047 error in position in error for L=2
%AttkScale=0.073; %i'm supposed to hit 0.08 error in position in error for L=3
AttkScale=25/0.7935; %i'm supposed to hit 0.106 error in position in error for L=4
%AttkScale=0.208; %i'm supposed to hit 0.179 error in position in error for L=6
refx = zeros(4,2000); % reference value for the controller



% [x_ea,y_ea,ya_ea,x_hat_ea,u_ea,a_ea] = SimulateClosedLoopSystemV2GAasEA_CDC18(A,B,C,Gatk,W,U,Q,R,x0,IsAttk,AtkDuration,AttkScale,refx,5);
% %beep
% 
% colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
% plot([0 200]*0.02,[0 0],'k--','LineWidth',2);
% hold on
% plot([0 200]*0.02,[0.5 0.5],'-.','LineWidth',2,'Color',colors(1,:));
% plot([0 200]*0.02,[0.3 0.3],'-.','LineWidth',2,'Color',colors(2,:));
% plot([1:200]*0.02,x_ea(1,:),'Color',colors(1,:))
% plot([1:200]*0.02,x_ea(3,:),'Color',colors(2,:))
% legend('reference','allowed distance error [m]','allowed angle error [rad]','Distance [m]','Angle of wheels [rad]','Orientation','vertical');
% grid
% xlabel('Time [s]');
% ylabel('Amplitude');
% ylim([-0.02 0.6]);
