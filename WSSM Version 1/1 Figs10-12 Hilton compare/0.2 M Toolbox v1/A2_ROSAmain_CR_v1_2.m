%% ROSA: Rate Of Spread Simulator - AAA method (Multiply Connected) -- CODE (DO NOT RUN)
function [bigZ, bigC, bigJ, merdata] = A2_ROSAmain_CR_v1_2(bigz,bigc,J,v0,delta,alpha,beta,lambda,tau,U,s,tstep,steps,resl,rkswit,pcswit,intswit,imswit)
%% Setup
lastwarn(''); warning('off','MATLAB:rankDeficientMatrix'); %reset old warnings
for j=1:J
    bigZ{j,1}=bigz{j}; bigC{j,1}=bigc{j}; bigJ{1}=J; %puts initial shapes into array bigZ
end
merdata = []; totaltime=0; 
 
%% The TimeStepping Loop
for m=1:steps
    tic;  % set timer

    % Fire Timestepping
    bigz = ROSAaaa_CR_v1_2(bigz,bigc,merdata,J,v0,delta,alpha,beta,lambda,tau,U,s,tstep,rkswit,pcswit,resl,intswit,imswit);

    % Fire Merger
    [bigz, bigc, merdata, J] = ROSAmerger_CR_v1_2(bigz,bigc,merdata,J,resl,intswit,imswit);

    % Update Lists
    for j=1:J
        bigZ{j,m+1}=bigz{j}; %update Z list
        bigC{j,m+1}=bigc{j}; %update C list
    end
    bigJ{m+1}=J;

    % Fire Image
    if imswit==1
        figure(1) %add to figure 1 on each timestep
        scsz = get(0,'ScreenSize'); set(gcf, 'Position',  [scsz(3)/20, scsz(4)/8, scsz(3)/2, scsz(3)/3])
        axl=4; axes('position',[.02 .3 .5 .5]);
        for k=1:m+1
            for j=1:bigJ{k}
                plot(real(bigZ{j,k}),imag(bigZ{j,k})), hold on,
            end
            if m==1 && k==1
                axis([-axl axl -axl axl]), daspect([1 1 1]), pause(0.05)
            end
        end
        hold off, axis([-axl axl -axl axl]), axis square, %daspect([1 1 1]),
        title("Wildfires",'interpreter','latex'),
    end

    % Update Time and User Message
    elapsedtime=toc; totaltime=totaltime+elapsedtime;
    if m~=1
          fprintf(repmat('\b',1,lineLength))
    end
    lineLength= fprintf('AAA-LS Fire: Step %d of %d completed in %.1f seconds. \n Total time: %.1f seconds.\n',m,steps,elapsedtime,totaltime); %feedback to user on time taken
end
end
