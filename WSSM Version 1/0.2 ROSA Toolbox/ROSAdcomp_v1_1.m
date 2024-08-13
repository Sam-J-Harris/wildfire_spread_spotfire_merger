function bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,prdt,shswt)
if shswt==2 || shswt==21
    load(prdt,'bigDataPack');
    bigZ0 = bigDataPack{1}; bigC0 = bigDataPack{2}; bigJ0 = bigDataPack{3}; merdata0 = bigDataPack{4}; tmax0 = bigDataPack{5};
    bigZ = [bigZ0 bigZ1(:,2:end)]; bigC = [bigC0 bigC1(:,2:end)]; bigJ = [bigJ0 bigJ1(2:end)];
    merdata = [merdata0, merdata1]; tmax = tmax0+tmax1;
    bigDataPack = {bigZ,bigC,bigJ,merdata,tmax};
else
    bigDataPack = {bigZ1,bigC1,bigJ1,merdata1,tmax1}; % big data pack for easy saving
end
end

% % Test
% prdt1 = 'bigDataPack_v1_2.mat'; prdt2 = 'bigDataPack_v1_2_2.mat';
% load(prdt1,'bigDataPack'); bigDataPacka = bigDataPack;
% load(prdt2,'bigDataPack'); bigDataPackb = bigDataPack;
% 
% bigZ0 = bigDataPacka{1}; bigC0 = bigDataPacka{2}; bigJ0 = bigDataPacka{3}; merdata0 = bigDataPacka{4};
% bigZ1 = bigDataPackb{1}; bigC1 = bigDataPackb{2}; bigJ1 = bigDataPackb{3}; merdata1 = bigDataPackb{4};
% bigZ = [bigZ0 bigZ1(:,2:end)]; bigC = [bigC0 bigC1(:,2:end)]; bigJ = [bigJ0 bigJ1(2:end)];
% merdata = [merdata0, merdata1];
% bigDataPack = {bigZ,bigC,bigJ,merdata};
