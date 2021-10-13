clc;
clear;
close all;

%load in joint configuration
q = importdata('q.csv').';
%visualize the configuration
robot = importrobot('robot.urdf');
robot.DataFormat = 'column';
show(robot,q,'visuals','on')