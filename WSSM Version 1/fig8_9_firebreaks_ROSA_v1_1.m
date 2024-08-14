%% ROSA: Rate Of Spread Simulator - AAA-LS method (Multiple Spotfires) -- RUN
% Figures 8 and 9 - firebreak examples
%
% Author's runtime: 
%   Fig8 (fgswt=8): 694.1 seconds.
%   Fig9a (fgswt=91): 110.5 seconds.
%   Fig9b (fgswt=92): 457.5 seconds.
%
% Code:
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off'); % clear workspace
fprintf('ROSA Activated: Beginning Spotfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
v0 = 1; alpha=0.5; delta=0.0; % basic ROS; rad/conv ratio (list of values); curvature param.
beta = 10; lambda = 2*beta; % pyrogenic and ambi wind params.
Umag=0.0; Uang=0; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % ambiwind magnitude and angle.

fgswt = 92; % figure choice: 8 = fig8, 91 = fig9a; 92 = fig9b
if fgswt==8 % three fires and a road
    alpha=0.25; % different alpha used 
    tstep=0.05; steps=72; spc=8;
    shswt = 124; rkswt=0; % Runge-Kutta (RK) switch = standard (0), RK2 (2) or RK4 (1 or 4) timestepping.
elseif fgswt==91 % circular fire and gap
    beta = 20; % different beta used
    tstep=0.1; steps=30; spc=3;
    shswt = 122; rkswt=0;
else % circular fire and lake
    tstep = 0.1; steps=40; spc=2;
    shswt = 125; rkswt=4;
end
 
tmin=0; tmaxe = tstep*steps; % min time, expected maximum time.
prdt = 'bigDataPack_example.mat'; % previous wildfire data (only for shswt=2).; % previous fire line data (only for shswt = 2 or 21, put 0 otherwise.
inswt=1; % interpolate polygon switch = off (0), on at each time step (1).
pcswt=1; % pole control switch = off (0), on (1).
imswt=0; % image display switch = off (0), on (1).

[resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt); % initial fire lines information - see function.

%% MAIN CODE AND PLOTTING
[bigZ1, bigC1, bigJ1, merdata1, tmax1, rtot1] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt,shswt);
bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt); % compile prdt and crdt into big data pack
%load('bigDataPack_example.mat'); shswt = 224; imswt=0; spc = 8; % bring in prev data (comment out if needed)
bigZ = bigDataPack{1}; bigJ = bigDataPack{3}; 
ROSAplot_v1_1(bigZ,bigJ,spc,shswt,imswt,1);
