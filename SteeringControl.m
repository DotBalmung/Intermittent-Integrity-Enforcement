% Singer model
clc;
clear all;
close all;

m = 1573; %kg
Iz = 2873; %vehicle yaw moment of inertia
lf = 1.1; %length in front of c.g.
lr = 1.58; % length behind of c.g. (rear)
Caf = 80000; %Cornering coefficient front
Car = 80000; %Cornering coefficient rear
Vx = 35; %mps
Ac = [0 1 0 0;0 -2*(Caf+Car)/(m*Vx) 2*(Caf+Car)/m 2*(-Caf*lf+Car*lr)/(m*Vx); 0 0 0 1 ; 0 -2*(Caf*lf-Car*lr)/(Iz*Vx) 2*(Caf*lf-Car*lr)/Iz -2*(Caf*lf^2+Car*lr^2)/(Iz*Vx)];
Bc = [0;2*Caf/m;0;2*Caf*lf/Iz];
Cc = eye(4);
Dc = 0;

Q = 0.001*eye(4);
R = 0.001*eye(4);
addProb = 0.0001;
FalsePosProb = 0.05;

%% Discretize the model with sampling rate Ts
Ts = 0.02;
csys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(csys,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%% plot the curve
% e1 = PlotErrorCurve(A,B,C,D,[],0,[],Q,R,FalsePosProb,[],[90 0],addProb,[]);
% hold on;
% e2 = PlotErrorCurve(A,B,C,D,[],0,[],Q,R,FalsePosProb,[],[90 5],addProb,[]);
% xlabel('sample');
% ylabel('Error norm ||e_k||_2');
% xlim([1 30]);
% grid on;
% legend('no policy','L=4');

% figure;
% e3 = PlotErrorCurve(A,B,C,D,[],1,[],Q,R,FalsePosProb,[],15,addProb,[]);
for i=1:5
    e3 = PlotErrorCurveWithF(A,B,C,D,[],1,[],Q,R,FalsePosProb,[],15,addProb,[],i);
    hold on;
end
legend('f=1','f=2','f=3','f=4','f=5');
% xlabel('Inter-enforcement distance l_3');
% ylabel('Maximal introduced error e_3^{max}');
% xlim([1 10]);
% seg1 = 1:10;
% point1 = [1 0];
% point2 = [10 e(10)];
% hold on;
% plotpoints = [point1;point2];
% plot(plotpoints(:,1),plotpoints(:,2),'--');
% legend('Original curve','linear approximation');
% %% compute slopes and offsets
% slope1 = (point2(2)-point1(2))/(point2(1)-point1(1))
% offset1 = point1(2)-slope1*point1(1)
