%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
% Figure 7 - three circular wildfires merging
%
% Author's runtime: 1656.5 seconds = 27.6 minutes.
%
% Parameter values:
%   v0 = 1; alpha=0.5; delta=0.0; 
%   beta = 20; lambda = 2*beta; 
%   Umag=0.0; Uang=pi/2; 
%   tstep=0.02; steps = 100; spc=10; (tstepm = 0.0005 - in ROSAtstep).
%   shswt=07; rkswt=4; inswt=1; pcswt=1; imswt=0; 
%
% Code:
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off'); % clear workspace
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alpha=0.5; delta=0.0; % basic ROS; rad/conv ratio; curvature param.
beta = 20; lambda = 2*beta; % pyrogenic and ambi wind params.
Umag=0.0; Uang=pi/2; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.
 
tstep=0.02; steps = 100; spc=10; % size of each step, number of steps, plot spacing (1: all timesteps plotted).
tmin=0; tmaxe = tstep*steps; % min and expected max time - may be different if emergency RK1 used.

shswt=07; % shape switch = circles (0), Hilton 2018 Fig 8 (1) - see function for more.
prdt = 'bigDataPack_example.mat'; % previous wildfire data (only for shswt=2+).
rkswt=4; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1) - pad polygon at each time step so that the resolution of points is fixed.
pcswt=1; % pole control switch = off (0), on (1) - manually remove Froissart doublets in AAA step.
imswt=0; % image display switch = off (0), on (1) - see images during the timestepping procedure.

[resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt); % initial fire line shapes - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1, tmax1, rtot1] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt,shswt); % main time stepping algorithm - see function.
bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt); % compile previous data (if applicable) and current data into big data pack - see function.
%load('bigDataPack_example.mat'); spc = 10; shswt = 21; imswt = 0; %bring in prev data (comment out)
bigZ = bigDataPack{1}; bigJ = bigDataPack{3}; % extract Z (boundary data) and J (no. of wildfires) values.
ROSAplot_v1_1(bigZ,bigJ,spc,shswt,imswt,1); % plot fire line evolution - see function.
