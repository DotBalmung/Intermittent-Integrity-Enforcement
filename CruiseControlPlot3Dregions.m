clc;clear all;close all;

% Load all of the system parameters. The file contains (among other info):
% System dynamics (A,B,C,D)
% Noise on states and the sensors (Q,R)
load('CruiseControlGA.mat');
%define colors for figures, so that 3D figure can be more legible
colors = ['y','r','b','k','m','g','c'];

%Find error regions for the system without integrity enforcement
[regs,norms] = FindErrorRegionsOurs(A,B,C,D,[1;0;0],[0;0;0],Q,R,7.82,0.001,zeros(1,4),0.05,5);
%Find error regions for the system with integrity enforcement on k=4
[regs2,norms2] = FindErrorRegionsOurs(A,B,C,D,[1;0;0],[0;0;0],Q,R,7.82,0.001,[0 0 0 1],0.05,5);

%Plot all of the figures
for i=1:4
Ellipse_plot(regs{i},zeros(length(A),1),colors(1+rem(i,7)));
hold on;
end
Ellipse_plot(regs2{4},zeros(length(A),1),colors(1+rem(5,7)));