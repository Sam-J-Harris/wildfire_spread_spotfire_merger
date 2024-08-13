%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
% Figure 6 - Two circular fires with various values of alpha and Umag
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off');
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic;

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alphaL = [0 0.5 1]; delta=0.0; % basic ROS; rad/conv ratio (list of values); curvature param.
beta = 1.5; lambda = 3; % pyrogenic and ambi wind params.
UmagL = [0 1]; Uang=pi/2; % ambiwind magnitude (list of values) and angle.

fgswt = 11; % choice of alpha and U -- 1: alpha=0, 2: alpha = 0.5, 3: alpha=1 -- 10: U=0, 20: U=1. 
alpha = alphaL(mod(fgswt,10)); Umag = UmagL(floor(fgswt/10)); U=Umag*cos(Uang)+1i*Umag*sin(Uang);

if floor(fgswt/10) == 2 % smaller timestep for Awind included
    tstep=0.005; steps = 50; spc=10;
else
    tstep=0.05; steps = 5; spc=1; % size of each step, number of steps, plot spacing (1: all timesteps plotted).
end
tmin=0; tmax = tstep*steps; tvec = linspace(tmin,tmax,steps+1); % min and max times, time vector.

shswt=06; % shape switch = circles (0), Hilton 2018 Fig 8 (1) - see function for more.
prdt = 'bigDataPack_v1_7_t1070.mat'; % previous wildfire data (only for shswt=2).; % previous fire line data (only for shswt = 2 or 21, put 0 otherwise.
rkswt=0; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1).
pcswt=1; % pole control switch = off (0), on (1).
imswt=0; % image display switch = off (0), on (1).

[resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt); % initial fire lines information - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1, tmax1] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt);
bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,prdt,shswt); % compile prdt and crdt into big data pack
%load('bigDataPack_v1_5_t160.mat'); spc = 1; shswt = 21; imswt = 0; %bring in prev data (comment out)
bigZ = bigDataPack{1}; bigJ = bigDataPack{3};
ROSAplot_v1_1(bigZ,bigJ,spc,shswt,imswt,1);