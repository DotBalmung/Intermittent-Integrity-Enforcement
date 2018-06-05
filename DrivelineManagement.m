% Singer model
clc;
clear all;
close all;
%% Another set of parameters
J_1=0.0974;
J_2=0.022;
k=5.5;
c=0;
b_1=1.7;
b_2 = 0.66
l=0;
i=1;

%% System matrices
Ac = [0 1/i -1; -k/i/J_1 -(b_1+c/i^2)/J_1 c/i/J_1; k/J_2 c/i/J_2 -(c+b_2)/J_2];
Bc = [0; 1/J_1; 0];
Cc = [0 1 0; 0 0 1];
Dc = 0;

%Q = 0.001*eye(3);
%R = 0.001*eye(2);
%addProb = 0.0001;
%FalsePosProb = 0.05;

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
    e3 = PlotErrorCurveMoreF(A,B,C,D,[],1,[],[],[],[],[],15,[],[],i);
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
