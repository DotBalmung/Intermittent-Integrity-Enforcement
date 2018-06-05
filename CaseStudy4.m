% Cruise control case study
clc;
clear all;
close all;

% define parameters related to the model
% Model is defined in continuous state-space
tau = 2.5;%s - time gap distance
Tl = 0.45; %eq (1), related to generalized vehicle longitudinal dynamic
Kl = 1;
Ac = [0 1 -tau ; 0 0 -1 ; 0 0 -1/Tl];
Bc = [0 0 ; 0 1 ; Kl/Tl 0];
Cc = eye(3);
Dc = 0;

Q = 0.01*eye(3);
R = 0.01*eye(3);

%% Discretize the model with sampling rate Ts
Ts = 0.02;
csys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(csys,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%Increase in probability of the detection induced by the attacker
probIncrease = 0.01;

%% plot the curve
e = PlotErrorCurve(A,B,C,D,[],1,[],Q,R,[],[],30,probIncrease,[]);
xlim([1 30]);
xlabel('Inter-enforcement distance l_1');
ylabel('Maximal introduced error e_1^{max}');

%% plot linear parts
seg1 = 1:2;
seg2 = 2:3;
seg3 = 3:18;
seg4 = 18:30;

p1 = polyfit(seg1',e(seg1),1);
p2 = polyfit(seg2',e(seg2),1);
p3 = polyfit(seg3',e(seg3),1);
p4 = polyfit(seg4',e(seg4),1);

hold on;
linfun1 = polyval(p1,seg1);
linfun2 = polyval(p2,seg2);
linfun3 = polyval(p3,seg3);
linfun4 = polyval(p4,seg4);

%% find continuous function
point1 = [1 linfun1(1)];
point2 = [2 linfun1(2)];
point3 = [3 linfun3(1)];
point4 = [18 linfun4(1)];
point5 = [30 linfun4(13)];

plotpoints = [point1;point2;point3;point4;point5];
plot(plotpoints(:,1),plotpoints(:,2),'--');
legend('Original curve','linear approximation');

%% compute slopes and offsets
slope1 = (point2(2)-point1(2))/(point2(1)-point1(1))
offset1 = point1(2)-slope1*point1(1)

slope2 = (point3(2)-point2(2))/(point3(1)-point2(1))
offset2 = point2(2)-slope2*point2(1)

slope3 = (point4(2)-point3(2))/(point4(1)-point3(1))
offset3 = point3(2)-slope3*point3(1)

slope4 = (point5(2)-point4(2))/(point5(1)-point4(1))
offset4 = point4(2)-slope4*point4(1)