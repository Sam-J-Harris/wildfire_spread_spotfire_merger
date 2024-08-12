%% ROSS: Rate Of Spread Simulator - ODE15i method (Single Wildfire) -- RUN
% Example code for calculating wildfire evolution of a single fire line.
% Change the shswt and shinp variables to choose initial fire line shape.
% Associated functions can be found in 0.1 ROSS Toolbox.
%
% Code:
clear; clc; close all; set(0,'defaultTextInterpreter','latex'); warning('off'); % clear workspace
fprintf('ROSS Activated: Beginning Wildfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 =1; alpha=0.5; delta=0.00; % basic ROS; rad/conv ratio; curvature param.
beta=10; lambda=20; % pyrogenic and ambi wind params.
Umag=1.0; Uang=0; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.

N=32; % Laurent series truncation
shswt=0; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape).
shinp = 21; % shape input -- see 0.1 ROSS Toolbox.

tmin=0; tmax=0.5; % min and max time.
steps = 10; spc = 1; % timesteps, plot spacing (1 => all timesteps plotted).
tvec=linspace(tmin,tmax,steps+1); % time vector.

%% MAIN CODE AND PLOTTING
[Z1, RE1, mRE1] = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,U); % main ODE solver
ROSSplot_v1_1({Z1},spc,1) % plotting function
totaltime=round(toc,2); fprintf("Fire Complete. Total time = "+num2str(totaltime)+" seconds.\n"); % stop timer, output message for the user.