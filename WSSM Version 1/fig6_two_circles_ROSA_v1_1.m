%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
% Figure 6 - two circular fires with various values of alpha and Umag
%
% Author's runtime: 
%   Fig6a (U=0, alpha=0.0, fgswt=11): 18.7 seconds.
%   Fig6b (U=0, alpha=0.5, fgswt=12): 18.4 seconds.
%   Fig6c (U=0, alpha=1.0, fgswt=13): 20.4 seconds.
%   Fig6d (U=1, alpha=0.0, fgswt=21): 202.8 seconds.
%   Fig6e (U=1, alpha=0.5, fgswt=22): 205.9 seconds.
%   Fig6f (U=1, alpha=1.0, fgswt=23): 221.1 seconds.
%
% Parameter values:
%   v0=1; alpha = [0 0.5 1]; delta=0.0; 
%   beta=1.5; lambda=3;
%   Umag=[0 1]; Uang=pi/2; 
%   tstep = [0.005 0.05]; steps = [50 5]; spc = [10 1]; rkswt = 0;
%
% Code:
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off'); % clear workspace
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alphaL = [0 0.5 1]; delta=0.0; % basic ROS; rad/conv ratio (list of values); curvature param.
beta = 1.5; lambda = 3; % pyrogenic and ambi wind params.
UmagL = [0 1]; Uang=pi/2; % ambiwind magnitude (list of values) and angle.

fgswt = 23; % choice of alpha and U -- 1: alpha=0, 2: alpha = 0.5, 3: alpha=1 -- 10: U=0, 20: U=1. 
alpha = alphaL(mod(fgswt,10)); Umag = UmagL(floor(fgswt/10)); U=Umag*cos(Uang)+1i*Umag*sin(Uang);

if floor(fgswt/10) == 2 % smaller timestep for Awind included
    tstep=0.005; steps = 50; spc=10;
else
    tstep=0.05; steps = 5; spc=1; % size of each step, number of steps, plot spacing (1: all timesteps plotted).
end
tmin=0; tmaxe = tstep*steps; % min and expected max time - may be different if emergency RK1 used.

shswt=06; % shape switch = circles (0), Hilton 2018 Fig 8 (1) - see function for more.
prdt = 'bigDataPack_example.mat'; % previous wildfire data (only for shswt=2+).
rkswt=0; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1) - pad polygon at each time step so that the resolution of points is fixed.
pcswt=1; % pole control switch = off (0), on (1) - manually remove Froissart doublets in AAA step.
imswt=0; % image display switch = off (0), on (1) - see images during the timestepping procedure.

[resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt); % initial fire line shapes - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1, tmax1,rtot1] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt); % main time stepping algorithm - see function.
bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt); % compile previous data (if applicable) and current data into big data pack - see function.
%load('bigDataPack_example.mat'); spc = 1; shswt = 21; imswt = 0; % bring in prev data (comment out if needed)
bigZ = bigDataPack{1}; bigJ = bigDataPack{3}; % extract Z (boundary data) and J (no. of wildfires) values.
ROSAplot_v1_1(bigZ,bigJ,spc,shswt,imswt,1); % plot fire line evolution - see function.