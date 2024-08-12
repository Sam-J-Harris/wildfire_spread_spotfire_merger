%% ROSS: Rate Of Spread Simulator - ODE15i method (Single Wildfire) -- RUN
% Figure 3 - pyrogenic wind comparison on thin ellipse with beta = 0.05, 5 
% and 10.
%
% Author's runtime: 623.32 seconds
%
% Parameter values:
%   v0=1; alpha = 0.5; delta=0.00; 
%   beta=0.05, 5, 10; lambda=10;  
%   Umag=1.0; Uang = pi/2;
%
%   N=128; shswt = 0; shinp = 12;
%   tmin = 0; tmax = 0.5; steps=1000; spc = 101;
%
% Code:
clear; clc; close all; set(0,'defaultTextInterpreter','latex'); warning('off'); % clear workspace
fprintf('ROSS Activated: Beginning Wildfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 =1; alpha=0.5; delta=0.00; % basic ROS; rad/conv ratio; curvature param.
beta=0.05; lambda=10; % pyrogenic and ambi wind params.
Umag=1.0; Uang=pi/2; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.

N=128; % series truncation
shswt=0; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape).
shinp = 12; % shape input -- see Toolbox.

tmin=0; tmax=0.5; % min and max time.
steps = 1000; spc = 101; % timesteps, plot spacing (1 => all timesteps plotted).
tvec=linspace(tmin,tmax,steps+1); % time vector.

%% MAIN CODE AND PLOTTING
[Z1, RE1, mRE1] = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,0.01,alpha,beta,lambda,U); % beta = 0.05, O(10^{-2}) curvature included.
[Z2, RE2, mRE2] = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,delta,alpha,5,lambda,U); % beta = 5.
[Z3, RE3, mRE3] = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,delta,alpha,10,lambda,U); % beta = 10.

ROSSplotPY_v1_1({Z1,Z2,Z3},spc,1); % plot all three wildfire evolutions
totaltime=round(toc,2); fprintf("Fire Complete. Total time = "+num2str(totaltime)+" seconds.\n"); % stop timer, output message for the user.