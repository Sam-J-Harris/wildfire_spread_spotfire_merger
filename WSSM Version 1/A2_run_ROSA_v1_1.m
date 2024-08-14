%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
% Example code for calculating wildfire evolution of multiple fire lines.
% Change the shswt variable to choose initial fire line shape.
% Associated functions can be found in 0.2 ROSA Toolbox.
%
% Code:
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off'); % clear workspace
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alpha=0.5; delta=0.0; % basic ROS; rad/conv ratio; curvature param.
beta = 10; lambda = 20; % pyrogenic and ambi wind params.
Umag=0.0; Uang=pi/2; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.

tstep=0.01; steps = 5; spc=1; % size of each time step, number of time steps, plot spacing (1: all timesteps plotted).
tmin=0; tmax = tstep*steps; tvec = linspace(tmin,tmax,steps+1); % min and max times, time vector.

shswt=0; % shape switch = circles (0), Hilton 2018 Fig 8 (1) - see function for more.
prdt = 'bigDataPack_v1_3.mat'; % previous wildfire data (only for shswt=2+).
rkswt=0; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1) - pad polygon at each time step so that the resolution of points is fixed.
pcswt=1; % pole control switch = off (0), on (1) - manually remove Froissart doublets in AAA step.
imswt=1; % image display switch = off (0), on (1) - see images during the timestepping procedure.

[resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt); % initial fire line shapes - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1,tmax1,rtot1] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt); % main time stepping algorithm - see function.
bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt); % compile previous data (if applicable) and current data into big data pack - see function.
bigZ = bigDataPack{1}; bigJ = bigDataPack{3}; % extract Z (boundary data) and J (no. of wildfires) values.
ROSAplot_v1_1(bigZ,bigJ,spc,shswt,imswt,1); % plot fire line evolution - see function.