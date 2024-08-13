function [bigZ, bigC, bigJ, merdata, tmax] = ROSAmain_v1_1(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,tstep,steps,resl,rkswt,pcswt,inswt,imswt)
% ROSA Main Code = produces fire line evolution and merger (if applicable)
%
% OUTPUTS
%   bigZ    = cell array of evolution of J fire lines over time
%   bigC    = cell array of fire line centres over time
%   bigJ    = cell array of number of wildfires J over time (decreases if
%               merger occurs.
%   merdata = information about mergers: provides initial overlapping fire
%               lines and final combined fire line.
%   tmax    = final time (varied if emergency RK1 timestepping used).
%
% END OF DOCUMENTATION
%
% Code:
%% Setup
lastwarn(''); warning('off','MATLAB:rankDeficientMatrix'); % reset old warnings
for j=1:size(bigz,2) % accounts for previous wildfire data 
    bigZ{j,1}=bigz{j}; bigC{j,1}=bigc{j}; bigJ{1}=J; % puts initial information into array bigZ, bigC and bigJ arrays.
end
merdata = {}; mcnt=size(bigz,2)-J; tmax=0; totaltime=0; % empty array ready for merge data, merge counter (starts at 0), start time at 0, set total runtime initially to zero.
 
%% The TimeStepping Loop
for m=1:steps
    tic; % set timer.
    [bigz,tmax] = ROSAtstep_v1_1(bigz,bigc,tmax,merdata,mcnt,J,v0,delta,alpha,beta,lambda,U,tstep,rkswt,pcswt,resl,imswt); % fire time step - see function.
    bigz = ROSAsmooth_v1_1(bigz,mcnt,J,resl,inswt,imswt); % fire line smoothing - see function.
    [bigz, bigc, merdata, mcnt, J] = ROSAmerger_v1_1(bigz,bigc,merdata,mcnt,J,resl,inswt,imswt); % fire merge (and smooth if necessary) - see function.

    for j=1:J % update lists.
        bigZ{j,m+1}=bigz{j}; %update Z list.
        bigC{j,m+1}=bigc{j}; %update C list.
    end
    bigJ{m+1}=J; % update J list.
    ROSAfimage_v1_1(bigZ,bigJ,m,imswt); % shows image of evolution at each time step, if imswt is turned on.

    elapsedtime=toc; totaltime=totaltime+elapsedtime; % update total runtime
    if m~=1, fprintf(repmat('\b',1,lineLength)), end
    lineLength= fprintf('AAA-LS Fire: Step %d of %d completed in %.1f seconds. \n Total time: %.1f seconds.\n',m,steps,elapsedtime,totaltime); %feedback to user on time taken for this time step.
end
end
