%% COMPARISON OF HILTON ET AL. 2018 FIGS 5, 7 AND 8
% Code:
clear; clc; close all; set(0,'defaultTextInterpreter','latex'); warning('off');

%% Parameters for both experiments
% Free to vary
delta = 0.01; % either 0.01 or 0.001
lambda = 2; % [2,8]
A = 0.24; % this is the main free parameter to match. OG: 0.625
alpha = 0.5; % [0.3,0.7] - try to get 0.5 or 0.6
rswt = 2; % switch: 1 = ROSSO only, 2 = ROSA only, 3 = both

% Single
N = 128;

% Multiple
steps8 = 800; msteps8 = 800; % steps =< msteps = max steps for whole thing
shswit=1; % circles (0), Hilton 2018 Fig 8 (1), prev Hilton data (2)
prevdata = 'v1_2_t12_alpha0.0_t0-600_workspace'; % bring in data (only for shswit=2)

% Fixed (or try to fix)
v0 = 0.1; % this must be 0.1
beta = lambda/2; % by Beer 1991
Umag = 1.0; Uang=pi/2; U=Umag*cos(Uang)+1i*Umag*sin(Uang); % same wind
tau = 0; s = fireterrain_1_1(1); % terrain params, never used
tmin = 0; % always start from 0
tmax5 = 0.3/A; tmax7 = 0.195/A; % scaled maximum times - single
mtmax8 = 0.84/A; tmax8 = (mtmax8/msteps8)*steps8; % actual max time
tstep8 = tmax8/steps8; % only to check that tstep isn't too high

load('crad_hilton18_fig5b_v7.mat','crad'); scl5=crad;
load('crad_hilton18_fig7e_csi_v7.mat','crad'); scl7=crad;
load('crad_hilton18_fig8e_csi_v7_1.mat','crad'); scl8=crad; % scales

%% Single Fires - Figs 5 and 7
if rswt == 1 || rswt == 3
fprintf('ROSSO_S Activated: Beginning Wildfire Simulation\n'); tic;
tvec5=linspace(tmin,tmax5,4); tvec7=linspace(tmin,tmax7,4); %time vectors

% Main Code
[Z1, RE1, mRE1] = ROSSOmain_v7_1(N,1,185,tvec5,v0,delta,alpha,beta/scl5,lambda/scl5,tau,U,s);
[Z2, RE2, mRE2] = ROSSOmain_v7_1(N,1,187,tvec7,v0,delta,alpha,beta/scl7,lambda/scl7,tau,U,s);

ROSSOplot_v7_1({Z1,Z2},1,1), ROSSOplotHI_v7_1(Z1,185,2), ROSSOplotHI_v7_1(Z2,187,3),
Stotaltime=round(toc,2); fprintf("Wildfire Complete. Total time = "+num2str(Stotaltime)+" seconds.\n");
end

%% Multiple Fires - Fig 8
if rswt==2 || rswt==3
fprintf('ROSA_M Activated: Beginning Spotfire Simulation\n');
tvec8 = linspace(tmin,tmax8,steps8+1); % time vector
imswit=0; % display images during the run
[resl, bigz, bigc, J] = ROSAshape_CR_v7_1(shswit,prevdata); % shape data

[bigZ1, bigC1, bigJ1, merdata] = A2_ROSAmain_CR_v1_2(bigz,bigc,J,v0,delta,alpha,beta/scl8,lambda/scl8,tau,U,s,tstep8,steps8,resl,0,1,1,imswit);
% for m=1:steps+1 % THIS IS IF J=1 INITIALLY IE FIRES HAVE MERGED IN THE PAST
%     bigZ1{2,m}=0; bigC1{2,m}=0;
% end
if shswit==2
    load(prevdata,'bigDataPack');
    bigZ0 = bigDataPack{1}; bigC0 = bigDataPack{2}; bigJ0 = bigDataPack{3}; merdata = bigDataPack{4};
    bigZ = [bigZ0 bigZ1(:,2:end)]; bigC = [bigC0 bigC1(:,2:end)]; bigJ = [bigJ0 bigJ1(2:end)];
    bigDataPack = {bigZ,bigC,bigJ,merdata};
    ROSAplot_CR_v7_1(bigZ,bigJ,100,shswit,4);
else
    bigDataPack = {bigZ1,bigC1,bigJ1,merdata};
    ROSAplot_CR_v7_1(bigZ1,bigJ1,100,shswit,4);
end
end