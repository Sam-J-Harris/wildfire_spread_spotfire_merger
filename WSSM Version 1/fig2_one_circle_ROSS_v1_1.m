%% ROSS: Rate Of Spread Simulator - ODE15i method (Single Wildfire) -- RUN
% Figure 2 - single circular wildfire in a uniform wind with relative error 
% plotted for N=32, 64, 128 Laurent series truncations.
%
% Author's runtime: 129.49 seconds.
%
% Parameter values:
%   v0=1; alpha = 0.5; delta=0.0; 
%   beta=10; lambda=20;
%   Umag=1.0; Uang=0; 
%   tmin = 0; tmax = 0.5; steps=1000; spc = 200
%   N = 32, 64 and 128
%
% Code:
clear; clc; close all; set(0,'defaultTextInterpreter','latex'); warning('off'); % clear workspace
fprintf('ROSS Activated: Beginning Wildfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 =1; alpha=0.5; delta=0.00; % basic ROS; rad/conv ratio; curvature param.
beta=10; lambda=20; % pyrogenic and ambi wind params.
Umag=1.0; Uang=0; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % wind properties

shswt=0; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape).
shinp = 21; % shape input -- see Toolbox.

tmin=0; tmax=0.5; % min and max time.
steps = 1000; spc = 200; % timesteps, plot spacing (1 => all timesteps plotted).
tvec=linspace(tmin,tmax,steps+1); % time vector.

%% MAIN CODE AND PLOTTING
[Z1, RE1, mRE1] = ROSSmain_v1_1(32,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,U); % N=32 run
[Z2, RE2, mRE2] = ROSSmain_v1_1(64,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,U); % N=64 run
[Z3, RE3, mRE3] = ROSSmain_v1_1(128,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,U); % N=128 run

bigDataPack = {Z3,tvec,RE1,RE2,RE3}; % put necessary final data into a "bigDataPack" - helps reduce storage space needed

% Z3 = bigDataPack{1}; tvec = bigDataPack{2}; spc = 200; % uncomment the next two lines if importing bigDataPack
% RE1 = bigDataPack{3}; RE2 = bigDataPack{4}; RE3 = bigDataPack{5}; 

ROSSplot_v1_1({Z3},spc,1) % plots the wildfire evolution
ROSSplotRE_v1_1(tvec,{RE1,RE2,RE3},2) % plots the comparison of the relative error for the N=32, 64 and 128 runs.

totaltime=round(toc,2); fprintf("Fire Complete. Total time = "+num2str(totaltime)+" seconds.\n"); % stop timer, output message for the user.