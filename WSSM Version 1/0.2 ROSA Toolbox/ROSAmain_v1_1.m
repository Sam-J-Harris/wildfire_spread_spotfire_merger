%% ROSA: Rate Of Spread Simulator - AAA method (Multiply Connected) -- CODE (DO NOT RUN)
function [bigZ, bigC, bigJ, merdata] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt)
% ROSA main code
% [bigZ, bigC, bigJ, merdata] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt)
%   = produces fire line evolution and merger (if applicable)
% 
% Inputs:
%
% Outputs:
%
% Code:
%% Setup
lastwarn(''); warning('off','MATLAB:rankDeficientMatrix'); % reset old warnings
for j=1:size(bigz,2) % accounts for prevdata 
    bigZ{j,1}=bigz{j}; bigC{j,1}=bigc{j}; bigJ{1}=J; % puts initial information into array bigZ, bigC and bigJ arrays.
end
merdata = {}; mcnt=0; totaltime=0; % empty array ready for merge data, merge counter (starts at 0), set totaltime initially to zero.
 
%% The TimeStepping Loop
for m=1:steps
    tic;  % set timer.
    bigz = ROSAtstep_v1_2(bigz,bigc,merdata,mcnt,J,v0,delta,alpha,beta,lambda,U,tstep,rkswt,pcswt,resl,imswt); % fire time step
    bigz = ROSAsmooth_v1_1(bigz,mcnt,J,resl,inswt,imswt); % fire line smoothing
    [bigz, bigc, merdata, mcnt, J] = ROSAmerger_v1_2(bigz,bigc,merdata,mcnt,J,resl,inswt,imswt); % fire merge (and smooth if nec)

    for j=1:J % update lists
        bigZ{j,m+1}=bigz{j}; %update Z list
        bigC{j,m+1}=bigc{j}; %update C list
    end
    bigJ{m+1}=J; % update J list
    ROSAfimage_v1_1(bigZ,bigJ,m,imswt); % shows image of evolution at each timestep, if selected

    elapsedtime=toc; totaltime=totaltime+elapsedtime; % update time
    if m~=1, fprintf(repmat('\b',1,lineLength)), end
    lineLength= fprintf('AAA-LS Fire: Step %d of %d completed in %.1f seconds. \n Total time: %.1f seconds.\n',m,steps,elapsedtime,totaltime); %feedback to user on time taken
end
end
