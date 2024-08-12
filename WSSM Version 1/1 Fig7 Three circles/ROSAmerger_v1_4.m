%% Appendix C: Fire Merger function
function [bigz, bigc, merdata, mcnt, J] = ROSAmerger_v1_4(bigz,bigc,merdata,mcnt,J,resl,inswt,imswt)
% = detects if fires are overlapping, merges any that do, smooths and interpolates each fire line.
% *some setup is needed in listing the fires: it will not detect whether e.g. fires 1 and 3 are overlapping
% Code:
imcnt = mcnt; % variable to assess if fires have merged this step
for j=1:J-1 % iterate through all fires
    if j<=J-1 % checks that J has not decreased during this loop
        mergeZ = firemerge(bigz{j}, bigz{j+1},bigz,j,J); % checks if fires j and j+1 are overlapping*
        if size(mergeZ,2)==1
            premerge = bigz; % original z data
            bigz{j} = mergeZ{1}; % replace fire line in jth array (centre of merged fire line = centre of jth fire line)
            if mcnt==1
                bigc{j} = 0.5-0.2i; % (centre of merged fire line = centre of second fire)
            end

            if j~=J-1
                for k = j+1:J-1 % shift each following fire line down by one entry
                    bigz{k} = bigz{k+1};
                    bigc{k} = bigc{k+1};
                end
            end
            J=J-1; mcnt=mcnt+1; % reduce J, add to the count of merged fires
            
            postmerge = bigz; % merged z data
            merdata{mcnt} = {premerge, postmerge}; % add to merge data (to see what happened on each merge step)
        else
            bigz{j}=mergeZ{1}; bigz{j+1}=mergeZ{2};
        end
    end
end

if mcnt>imcnt % if fires have merged this step
    bigz = ROSAsmooth_v1_4(bigz,mcnt,J,resl,inswt,imswt);
end
end

%% Appendix C1: Merge Function
function mergeZ = firemerge(z1,z2,bigz,j,J)
% = determine if fire lines z1 and z2 overlap and merge them with union ftn if they do.
% Code:
x1 = real(z1); y1 = imag(z1); x1 = x1(~isnan(x1)); y1 = y1(~isnan(y1));
x2 = real(z2); y2 = imag(z2); x2 = x2(~isnan(x2)); y2 = y2(~isnan(y2)); % remove any NaN points from data

poly1 = polyshape(x1,y1); poly2 = polyshape(x2,y2); % convert to polyshape objects

if overlaps(poly1,poly2) % determine if the fire lines overlap
    polyun = union(poly1,poly2); [X,Y] = boundary(polyun); Z = X+1i*Y; % create union (merge) of the two lines and extract the boundary
    if ~ispolycw(X,Y), Z = flip(Z); end % flip if anticlockwise
    Pts = size(Z,1); smno=20; Zsm = firesmooth(Z,Pts,smno); % smooth data

    % figure(6)
    % plot(real(Z),imag(Z)); daspect([1 1 1])

    curvZ = LineCurv4(Zsm); cmax = max(abs(curvZ)); pos = curvZ==cmax; % finds curvature of each point and identifies the maximum
    zmpk = Zsm(pos); zdis = abs(Zsm-zmpk.*ones(size(Zsm))); % identifies the peak of the merge: "zmpeak"
    if J~=2
        z3 = bigz{j+2}; zdis2 = abs(Zsm-z3(1)*ones(size(Zsm)));
        zdis = zdis+zdis2;
    end

    zdmax = max(zdis); % finds distance of each other pt from the merge peak
    
    while zdis(1)<zdmax % shifts Z coordinates until the start pt is the pt furthest from the merge
        Zsm = circshift(Zsm,-1); zdis = circshift(zdis,-1);
    end 
    mergeZ{1} = Zsm; 
else
    mergeZ{1} = z1; mergeZ{2} = z2; % output the input data if no overlap
end
end

%% Appendix C2: Fire Line Smooth and Interp Function
function znew = firesmooth(zold, Pts, smno)
% = smooths fire line and interpolates to desired numbers of pts
% Code:
znew = zold;
for k = 1:smno % smooths a certain number of times
    znew = smoothdata(znew,'gaussian', 4); znew = [znew; znew(1)]; % smooth and close
    znew = fireinterp(znew,Pts,'spline'); % interpolate smoothed data
end
end