function [resl, bigz, bigc, J] = ROSAshape_CR_v1_2(shswit,prevdata)
    resl = 100; 
    if shswit ==0
        sepd = 1.5; % resolution = pts/length, initial separation of fires
        A1 = 1; B1 = 0.0; c1 = -sepd; nsym1 = 2; % radius and symmetry of fire 1
        A2 = 1; B2 = 0.0; c2 = sepd; nsym2 = 2; % radius and symmetry of fire 2

        A = [A1, A2]; B = [B1, B2]; nsym = [nsym1, nsym2]; %big lists
        bigz = {}; bigc = {c1, c2}; J=size(bigc,2); % big arrays, J = no. of disks

        for j=1:J
            Pts = round(2*pi*resl*A(j)); tht = linspace(0,2*pi,Pts); zta = exp(1i*tht);
            bigz{j} = ((A(j)*zta.^(-1)+B(j)*zta.^(nsym(j)-1)+bigc{j}).');
        end
    elseif shswit==1
        load('bigz_hilton2018_fig8_csi_v6.mat', 'bigz');
        load('bigc_hilton2018_fig8_csi_v6.mat', 'bigc'); J = size(bigc,2); 
    else
        load(prevdata,'bigDataPack');
        Jtemp = bigDataPack{3}; J = Jtemp{end};
        Ztemp = bigDataPack{1}; Ctemp = bigDataPack{2};
        for j=1:J
             bigz{j} = Ztemp{j,end}; bigc{j} = Ctemp{j,end};
        end
    end
end