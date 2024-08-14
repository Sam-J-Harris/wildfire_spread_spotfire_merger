%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off');
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic;

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alpha=0.5; delta=0.0; % basic ROS; rad/conv ratio; curvature param.
beta = 20; lambda = 2*beta; % pyrogenic and ambi wind params.
Umag=0.0; Uang=pi/2; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.
 
tstep=0.02; steps = 100; spc=10; % size of each step, number of steps, plot spacing (1: all timesteps plotted).
tmin=0; tmaxe = tstep*steps; % min time, expected maximum time.

shswt=07; % shape switch = circles (0), Hilton 2018 Fig 8 (1) - see function for more.
prdt = 'fix1_v4_t200.mat'; % previous wildfire data (only for shswt=2).; % previous fire line data (only for shswt = 2 or 21, put 0 otherwise.
rkswt=4; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1).
pcswt=1; % pole control switch = off (0), on (1).
imswt=0; % image display switch = off (0), on (1).

[resl, bigz, bigc, J] = ROSAshape_v1_4(shswt,prdt); % initial fire lines information - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1, tmax1, rtot1] = ROSAmain_v1_4(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt);
bigDataPack = ROSAdcomp_v1_4(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt); % compile previous data (if applicable) and current data into big data pack.
%load('fix1_v6_t100.mat'); spc = 10; shswt = 21; imswt = 0; %bring in prev data (comment out)
bigZ = bigDataPack{1}; bigJ = bigDataPack{3};
ROSAplot_v1_4(bigZ,bigJ,spc,shswt,imswt,1);
