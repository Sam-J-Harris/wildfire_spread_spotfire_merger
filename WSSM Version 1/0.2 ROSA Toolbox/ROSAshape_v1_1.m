function [resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prdt)
% ROSA shape function
% [resl, bigz, bigc, J] = ROSAshape_v1_1(shswt,prevdata)
%   = generates initial fire line data.
% 
% Inputs:
%   shswt = shape switch: circles (0), Hilton 2018 Fig 8 (1), prev Hilton data (2).
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
    sepd = 1.1; % initial separation of fires
    A1 = 1; B1 = 0.0; c1 = -sepd; nsym1 = 2; % radius and symmetry of fire 1
    A2 = 1; B2 = 0.0; c2 = sepd; nsym2 = 2; % radius and symmetry of fire 2

    A = [A1, A2]; B = [B1, B2]; nsym = [nsym1, nsym2]; % big lists.
    bigz = {}; bigc = {c1, c2}; J=size(bigc,2); % set up z and c arrays, J = no. of fires.

elseif shswt == 07 % four fires, no merger (figure 7).

elseif shswt == 08 % three circular fires (figure 8).
    A1 = 0.5; B1 = 0.0; c1 = -1; nsym1 = 2; % radius and symmetry of fire 1
    A2 = 0.5; B2 = 0.0; c2 = 0.5-0.2i; nsym2 = 2; % radius and symmetry of fire 2
    A3 = 0.5; B3 = 0.0; c3 = 2.25+0.3i; nsym3 = 2; % radius and symmetry of fire 2

    A = [A1, A2, A3]; B = [B1, B2, B3]; nsym = [nsym1, nsym2, nsym3]; % big lists.
    bigz = {}; bigc = {c1, c2, c3}; J=size(bigc,2); % set up z and c arrays, J = no. of fires.

else % default to two circular fires otherwise
    sepd = 1.1; % initial separation of fires
    A1 = 1; B1 = 0.0; c1 = -sepd; nsym1 = 2; % radius and symmetry of fire 1
    A2 = 1; B2 = 0.0; c2 = sepd; nsym2 = 2; % radius and symmetry of fire 2

    A = [A1, A2]; B = [B1, B2]; nsym = [nsym1, nsym2]; % big lists.
    bigz = {}; bigc = {c1, c2}; J=size(bigc,2); % set up z and c arrays, J = no. of fires.
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
    for j=1:J
        Pts = round(2*pi*resl*A(j)); tht = linspace(0,2*pi,Pts); zta = exp(1i*tht);
        bigz{j} = ((A(j)*zta.^(-1)+B(j)*zta.^(nsym(j)-1)+bigc{j}).'); % update z array.
    end
end
end