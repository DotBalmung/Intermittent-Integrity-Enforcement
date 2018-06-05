% Differential braking
clc;
clear all;
close all;

m = 1820; %kg
Iz = 2922; %vehicle yaw moment of inertia
u = 27.78; % forward speed in mps
a = 1.07; %length in front of c.g.
b = 1.62; % length behind of c.g. (rear)
Ca1 = -141670; %Cornering coefficient front N/rad
Ca2 = -119700; %Cornering coefficient rear
T = 1.52; %m - track width
tau = 0.1; % hydraulic response time constant - right above (13)
Ac = [1/tau 0 0 0 0 ; 0 (Ca1+Ca2)/(m*u) (a*Ca1-b*Ca2-m*u^2)/(m*u) 0 0 ; T/(2*Iz) (a*Ca1-b*Ca2)/(Iz*u) (a^2*Ca1+b^2*Ca2)/(Iz*u) 0 0 ; 0 (Ca1+Ca2)/(m*u) (a*Ca1-b*Ca2)/(m*u) 0 0 ; 0 0 0 1 0];
Bc = [1/tau ; 0 ; 0 ; 0 ; 0];
Cc = [0 0 0 0 1];
Dc = 0;

Q = 0.01*eye(5);
R = 0.01;

%% Discretize the model with sampling rate Ts
Ts = 0.02;
csys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(csys,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%% plot the curve
e = PlotErrorCurve(A,B,C,D,[],1,[],[],[],[],[],10,[],[]);
xlim([1 10]);
xlabel('Inter-enforcement distance l_2');
ylabel('Maximal introduced error e_2^{max}');

%% plot linear parts
% polifyt didn't provide satisfiable results, so linear segments were
% chosen manualy
%
% seg1 = 1:2;
% seg2 = 2:10;
% 
% p1 = polyfit(seg1',e(seg1),1);
% p2 = polyfit(seg2',e(seg2),1);
% 
hold on;
% linfun1 = polyval(p1,seg1);
% linfun2 = polyval(p2,seg2);

%% find continuous function
% point1 = [1 linfun1(1)];
% point2 = [2 linfun2(1)];
% point3 = [10 linfun2(8)];
point1 = [1 0];
point2 = [2 219645];
point3 = [10 334468];

plotpoints = [point1;point2;point3];
plot(plotpoints(:,1),plotpoints(:,2),'--');
legend('Original curve','linear approximation');

%% compute slopes and offsets
slope1 = (point2(2)-point1(2))/(point2(1)-point1(1))
offset1 = point1(2)-slope1*point1(1)

slope2 = (point3(2)-point2(2))/(point3(1)-point2(1))
offset2 = point2(2)-slope2*point2(1)