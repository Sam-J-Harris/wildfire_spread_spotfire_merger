%% ROSS: Rate Of Spread Simulator - ODE15i method (Single Wildfire) -- RUN
% Figure 4 - evolution of various Laurent and non-Laurent shapes
%
% Author's runtimes: 
%   irregular pentagon:
%   hourglass:
%   bean:
%   blade:
%
% Code:
clear; clc; close all; set(0,'defaultTextInterpreter','latex'); warning('off'); % clear workspace
fprintf('ROSS Activated: Beginning Wildfire Simulation\n'); tic; % startup message to the user

%% PARAMETER SET-UP -- USER INPUT
% Common Parameters for all experiments
v0=1; alpha=0.5; Umag=1.0; % basic ROS; rad/conv ratio; ambient wind magnitude.
tmin=0; steps=1000; spc=200; N=128; % minimum time; number of timesteps computed; spacing between plotted isochrones; series truncation.

% Irregular pentagon
delta1=0.01; beta1=10; lambda1=20; % curvature; pyrowind; ambiwind;
Uang1=pi/6; U1 = Umag*cos(Uang1)+1i*Umag*sin(Uang1); % ambiwind direction.
shswt1=0; shinp1=51; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape); shape input -- see Toolbox.
tmax1=0.25; tvec1=linspace(tmin,tmax1,steps+1); % maximum time; time vector.

% Hourglass
delta2=0.05; beta2=18; lambda2=25; % curvature; pyrowind; ambiwind;
Uang2=-pi/4; U2=Umag*cos(Uang2)+1i*Umag*sin(Uang2); % ambiwind direction.
shswt2=0; shinp2 = 71; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape); shape input -- see Toolbox.
tmax2=0.4; tvec2=linspace(tmin,tmax2,steps+1); % maximum time; time vector.

% Bean
delta3=0.1; beta3=10; lambda3=20; % curvature; pyrowind; ambiwind;
Uang3=5*pi/6; U3=Umag*cos(Uang3)+1i*Umag*sin(Uang3); % ambiwind direction.
shswt3=1; shinp3 = 101; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape); shape input -- see Toolbox.
tmax3=0.2; tvec3=linspace(tmin,tmax3,steps+1); % maximum time; time vector.

% Blade
delta4=0.15; beta4=6; lambda4=15; % curvature; pyrowind; ambiwind;
Uang4=pi/3; U4=Umag*cos(Uang4)+1i*Umag*sin(Uang4); % ambiwind direction.
shswt4=1; shinp4 = 102; % shape switch = 0 (Laurent shape) , 1 (non-Laurent shape); shape input -- see Toolbox.
tmax4=0.2; tvec4=linspace(tmin,tmax4,steps+1); % maximum time; time vector.

%% MAIN CODE AND PLOTTING
% Main ODE solver and plotting - uncomment desired example.
%[Z1, RE1, mRE1] = ROSSmain_v1_1(N,shswt1,shinp1,tvec1,v0,delta1,alpha,beta1,lambda1,U1); ROSSplot_v1_1({Z1},spc,1); % irregular pentagon.
%[Z2, RE2, mRE2] = ROSSmain_v1_1(N,shswt2,shinp2,tvec2,v0,delta2,alpha,beta2,lambda2,U2); ROSSplot_v1_1({Z2},spc,2); % hourglass.
%[Z3, RE3, mRE3] = ROSSmain_v1_1(N,shswt3,shinp3,tvec3,v0,delta3,alpha,beta3,lambda3,U3); ROSSplot_v1_1({Z3},spc,3); % bean.
[Z4, RE4, mRE4] = ROSSmain_v1_1(N,shswt4,shinp4,tvec4,v0,delta4,alpha,beta4,lambda4,U4); ROSSplot_v1_1({Z4},spc,4); % blade.
totaltime=round(toc,2); fprintf("Fire Complete. Total time = "+num2str(totaltime)+" seconds.\n"); % stop timer, output message for the user.