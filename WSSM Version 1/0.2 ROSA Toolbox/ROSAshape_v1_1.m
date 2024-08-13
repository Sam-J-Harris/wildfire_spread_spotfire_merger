function [resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt)
% ROSA shape function
% [resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prevdata)
%   = generates initial fire line data.
% 
% Inputs:
%   shswt = shape switch: circles (0), Hilton 2018 Fig 8 (1), prev Hilton data (2),
%                           figure 6 two circles (06), figure 7 four fires (07), 
%                           figure 8 three circles (08).
%
% Outputs:
%   resl = resolution of all fire lines (points/length).
%   bigz = cell array of size J of initial fire line data (z coordinates) for the J fires.
%   bigc = cell array of size J of initial fire line centres for the J fires.
%   J = number of (initial) wildfires.
%
% Code:
resl = 100; % set resl = 100 for all experiments (can be changed).

if shswt == 06 % two circular fires for alpha comparison (figure 6).
    sepd = 1.35; % initial separation of fires
    A{1} = 1; B{1} = 0.0; bigc{1} = -sepd; nsym{1} = 2; % radius and symmetry of fire 1
    A{2} = 1; B{2} = 0.0; bigc{2} = sepd; nsym{2} = 2; % radius and symmetry of fire 2

elseif shswt == 07 % four fires, no merger (figure 7).
    A{1} = 0.5; B{1} = 0.0; bigc{1} = -1.5; nsym{1} = 2; % radius and symmetry of fire 1
    A{2} = 1; B{2} = 0.1; bigc{2} = 0+2i; nsym{2} = 5; % radius and symmetry of fire 2
    A{3} = 0.75; B{3} = 0.15; bigc{3} = 1.5; nsym{3} = 4; % radius and symmetry of fire 3
    A{4} = 1; B{4} = 0.5; bigc{4} = -2i; nsym{4} = 2; % radius and symmetry of fire 4

elseif shswt == 08 % three circular fires (figure 8).
    A{1} = 0.5; B{1} = 0.0; bigc{1} = -1; nsym{1} = 2; % radius and symmetry of fire 1
    A{2} = 0.5; B{2} = 0.0; bigc{2} = 0.5-0.2i; nsym{2} = 2; % radius and symmetry of fire 2
    A{3} = 0.5; B{3} = 0.0; bigc{3} = 2.5+0.3i; nsym{3} = 2; % radius and symmetry of fire 2

elseif shswt == 085 % variation of three circles
    A{1} = 0.5; B{1} = 0.0; bigc{1} = -1; nsym{1} = 2; % radius and symmetry of fire 1
    A{2} = 0.5; B{2} = 0.0; bigc{2} = 0.5-0.2i; nsym{2} = 2; % radius and symmetry of fire 2
    A{3} = 0.5; B{3} = 0.0; bigc{3} = 3.0+0.3i; nsym{3} = 2; % radius and symmetry of fire 2

else % default to two circular fires otherwise
    sepd = 1.1; % initial separation of fires
    A{1} = 1; B{1} = 0.0; bigc{1} = -sepd; nsym{1} = 2; % radius and symmetry of fire 1
    A{2} = 1; B{2} = 0.0; bigc{2} = sepd; nsym{2} = 2; % radius and symmetry of fire 2
end

if shswt==1 % Hilton et al. 2018 figure 8 spotfires.
    load('bigz_hilton2018_fig8_csi_v6.mat', 'bigz');
    load('bigc_hilton2018_fig8_csi_v6.mat', 'bigc'); J = size(bigc,2); 

elseif shswt==2 % load in previous Hilton wildfire data.
    load(prdt,'bigDataPack');
    Jtemp = bigDataPack{3}; J = Jtemp{end}; % need "temp" files if fires merged previously.
    Ztemp = bigDataPack{1}; Ctemp = bigDataPack{2};
    for j=1:J
         bigz{j} = Ztemp{j,end}; bigc{j} = Ctemp{j,end};
    end
    J = Jtemp{end};

elseif shswt==21 % load in previous n-shape wildfire data.
    load(prdt,'bigDataPack');
    Jtemp = bigDataPack{3}; Jst = Jtemp{1}; J = Jtemp{end}; % need "temp" files if fires merged previously.
    Ztemp = bigDataPack{1}; Ctemp = bigDataPack{2};
    for j=1:Jst
         bigz{j} = Ztemp{j,end}; bigc{j} = Ctemp{j,end};
    end

else % n-fold shapes
    J=size(bigc,2); % set up z and c arrays, J = no. of fires.
    for j=1:J
        Pts = round(2*pi*resl*A{j}); tht = linspace(0,2*pi,Pts); zta = exp(1i*tht);
        bigz{j} = ((A{j}*zta.^(-1)+B{j}*zta.^(nsym{j}-1)+bigc{j}).'); % update z array.
    end
end
end