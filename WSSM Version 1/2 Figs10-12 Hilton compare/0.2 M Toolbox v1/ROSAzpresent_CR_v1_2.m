% Presenting Solution
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off');

prevdata = 'v1_2_t11_alpha0.0_t0-600_workspace.mat'; % bring in data
%prevdata = 'v1_2_v04_t10.5_alpha0.5_t0-800_workspace.mat'; % bring in data
load(prevdata,'bigDataPack');
bigZ1 = bigDataPack{1}; bigC1 = bigDataPack{2}; bigJ1 = bigDataPack{3};
spc = 100; shswit = 2; %spc=20;
ROSAplot_CR_v1_2(bigZ1,bigJ1,spc,shswit,1);