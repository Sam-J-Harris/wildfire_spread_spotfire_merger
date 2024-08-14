%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off');
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic;

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alpha=0.5; delta=0.0; % basic ROS; rad/conv ratio; curvature param.
beta = 10; lambda = 2*beta; % pyrogenic and ambi wind params.
Umag=0.0; Uang=0; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.
 
tstep=0.01; steps = 10; spc=1; % size of each step, number of steps, plot spacing (1: all timesteps plotted).
tmin=0; tmaxe = tstep*steps; % min time, expected maximum time.

shswt=125; % shape switch. Fig 8 (road) = 124; Fig9a (gap) = 122; Fig9b (lake) = 125
prdt = 'lake_t30_v1.mat'; % previous wildfire data (only for shswt=2).; % previous fire line data (only for shswt = 2 or 21, put 0 otherwise.
rkswt=4; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1).
pcswt=1; % pole control switch = off (0), on (1).
imswt=0; % image display switch = off (0), on (1).

[resl, bigz, bigc, J] = ROSAshape_fbk_v1_1(shswt,prdt); % initial fire lines information - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1, tmax1, rtot1] = ROSAmain_fbk_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt,shswt);
bigDataPack = ROSAdcomp_fbk_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt); % compile prdt and crdt into big data pack
%load('tcr_t72_v1.mat'); 
shswt = 224; imswt=0; spc = 8;
bigZ = bigDataPack{1}; bigJ = bigDataPack{3}; 
ROSAplot_fbk_v1_1(bigZ,bigJ,spc,shswt,imswt,1);
