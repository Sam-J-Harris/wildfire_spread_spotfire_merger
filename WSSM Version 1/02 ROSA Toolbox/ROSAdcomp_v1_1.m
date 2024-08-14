function bigDataPack = ROSAdcomp_v1_1(bigZ1, bigC1, bigJ1,merdata1,tmax1,rtot1,prdt,shswt)
% ROSA data compilation function
%   = compiles previous wildfire data (if applicable) and current wildfire
%   data into one bigDataPack for easy saving and later use
%
%   Contains bigZ, bigC, bigJ, merdata, tmax and rtot data.
%
% Code:
if shswt==2 || shswt==21 % adding previous data to the data pack (if applicable).
    load(prdt,'bigDataPack');
    bigZ0 = bigDataPack{1}; bigC0 = bigDataPack{2}; bigJ0 = bigDataPack{3}; 
    merdata0 = bigDataPack{4}; tmax0 = bigDataPack{5}; rtot0 = bigDataPack{6};
    bigZ = [bigZ0 bigZ1(:,2:end)]; bigC = [bigC0 bigC1(:,2:end)]; bigJ = [bigJ0 bigJ1(2:end)]; % joins previous and current data together.
    merdata = [merdata0, merdata1]; tmax = tmax0+tmax1; rtot = rtot0+rtot1;
    bigDataPack = {bigZ,bigC,bigJ,merdata,tmax,rtot};
else
    bigDataPack = {bigZ1,bigC1,bigJ1,merdata1,tmax1,rtot1}; % big data pack for easy saving
end
end