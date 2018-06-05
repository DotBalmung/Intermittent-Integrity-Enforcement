clear all;close all;clc;


%% plot all the figures that show effects of attack on QoC

%% Cruise control : L13 - L5 - Lvar
figure(1);
%
subplot(1,3,1);
load('CruiseControlLis13.mat')
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
plot([0 2000]*0.02,[0 0],'k--','LineWidth',2);
hold on
plot([0 2000]*0.02,[1 1],'-.','LineWidth',2,'Color',colors(1,:));
plot([0 2000]*0.02,[0.5 0.5],'-.','LineWidth',2,'Color',colors(2,:));
plot([1:2000]*0.02,x_ea(1,:),'Color',colors(1,:))
plot([1:2000]*0.02,x_ea(2,:),'Color',colors(2,:))
legend('reference','allowed mean postion error [m]','allowed mean speed error [m/s]','Error in position [m]','Error in speed [m/s]','Orientation','horizontal');
grid
xlabel('Time [s]');
ylabel('Amplitude');
ylim([-0.21 1.1]);
%yticks([0 0.25 0.5 0.75 1]);
%
clear all;
subplot(1,3,2);
load('CruiseControlLis5.mat')
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
plot([0 2000]*0.02,[0 0],'k--','LineWidth',2);
hold on
plot([0 2000]*0.02,[1 1],'-.','LineWidth',2,'Color',colors(1,:));
plot([0 2000]*0.02,[0.5 0.5],'-.','LineWidth',2,'Color',colors(2,:));
plot([1:2000]*0.02,x_ea(1,:),'Color',colors(1,:))
plot([1:2000]*0.02,x_ea(2,:),'Color',colors(2,:))
%legend('reference','allowed mean postion error [m]','allowed mean speed error [m/s]','Error in position [m]','Error in speed [m/s]');
grid
xlabel('Time [s]');
%ylabel('Amplitude');
ylim([-0.21 1.1]);
%yticks([0 0.25 0.5 0.75 1]);
%
clear all;
subplot(1,3,3);
load('CruiseControlChangingLs.mat');
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
plot([0 2000]*0.02,[0 0],'k--','LineWidth',2);
hold on
plot([0 2000]*0.02,[1 1],'-.','LineWidth',2,'Color',colors(1,:));
plot([0 2000]*0.02,[0.5 0.5],'-.','LineWidth',2,'Color',colors(2,:));
plot([1:2000]*0.02,x_ea(1,:),'Color',colors(1,:))
plot([1:2000]*0.02,x_ea(2,:),'Color',colors(2,:))
%legend('reference','allowed mean postion error [m]','allowed mean speed error [m/s]','Error in position [m]','Error in speed [m/s]');
grid
xlabel('Time [s]');
%ylabel('Amplitude');
ylim([-0.21 1.1]);
%yticks([0 0.25 0.5 0.75 1]);

%% Lane position controller : L6 - L4 - Lvar
figure(2);
%
subplot(1,3,1);
load('LanePositionLis6.mat')
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
plot([0 2000]*0.02,[0 0],'k--','LineWidth',2);
hold on
plot([0 2000]*0.02,[0.2 0.2],'-.','LineWidth',2,'Color',colors(1,:));
plot([0 2000]*0.02,[0.18 0.18],'-.','LineWidth',2,'Color',colors(2,:));
plot([1:2000]*0.02,x_ea(1,:),'Color',colors(1,:))
plot([1:2000]*0.02,x_ea(3,:),'Color',colors(2,:))
legend('reference','allowed mean distance error [m]','allowed mean angle error [rad]','Distance from road center [m]','Angle of wheels to the road [rad]','Orientation','horizontal');
grid
xlabel('Time [s]');
ylabel('Amplitude');
ylim([-0.02 0.21]);
%yticks([0 0.05 0.1 0.15 0.2]);
%
subplot(1,3,2);
load('LanePositionLis4.mat')
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
plot([0 2000]*0.02,[0 0],'k--','LineWidth',2);
hold on
plot([0 2000]*0.02,[0.2 0.2],'-.','LineWidth',2,'Color',colors(1,:));
plot([0 2000]*0.02,[0.18 0.18],'-.','LineWidth',2,'Color',colors(2,:));
plot([1:2000]*0.02,x_ea(1,:),'Color',colors(1,:))
plot([1:2000]*0.02,x_ea(3,:),'Color',colors(2,:))
%legend('reference','allowed mean distance error [m]','allowed mean angle error [rad]','Distance from road center [m]','Angle of wheels to the road [rad]','Orientation','horizontal');
grid
xlabel('Time [s]');
%ylabel('Amplitude');
ylim([-0.02 0.21]);
%yticks([0 0.05 0.1 0.15 0.2]);
%
subplot(1,3,3);
load('LanePositionChangingLs.mat')
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
plot([0 2000]*0.02,[0 0],'k--','LineWidth',2);
hold on
plot([0 2000]*0.02,[0.2 0.2],'-.','LineWidth',2,'Color',colors(1,:));
plot([0 2000]*0.02,[0.18 0.18],'-.','LineWidth',2,'Color',colors(2,:));
plot([1:2000]*0.02,x_ea(1,:),'Color',colors(1,:))
plot([1:2000]*0.02,x_ea(3,:),'Color',colors(2,:))
%legend('reference','allowed mean distance error [m]','allowed mean angle error [rad]','Distance from road center [m]','Angle of wheels to the road [rad]','Orientation','horizontal');
grid
xlabel('Time [s]');
%ylabel('Amplitude');
ylim([-0.02 0.21]);
%yticks([0 0.05 0.1 0.15 0.2]);
