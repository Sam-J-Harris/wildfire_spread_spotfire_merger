%% Appendix C: Fire Merger function
function [bigz, bigc, merdata, mcnt, J] = ROSAmerger_v1_1(bigz,bigc,merdata,mcnt,J,resl,inswt,imswt,shswt)
% = detects if fire lines are overlapping, merges any that do, smooths and interpolates each fire line.
% *some setup is needed in listing the fires: it will not detect whether e.g. fires 1 and 3 are overlapping
%
% Code:
imcnt = mcnt; % "initial merge counter": variable to assess if fires have merged this step
for j=1:J-1 % iterate through all fires
    if j<=J-1 % checks that J has not decreased during this loop
        mergeZ = firemerge(bigz{j}, bigz{j+1},bigz,j,J,shswt); % checks if fires j and j+1 are overlapping*
        if size(mergeZ,2)==1 % if a merge has happened
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
            J=J-1; mcnt=mcnt+1; % reduce J by 1, add 1 to the count of merged fires "mcnt"
            
            postmerge = bigz; % merged z data
            merdata{mcnt} = {premerge, postmerge}; % add to merge data (to see what happened on each merge step)
        else % if no merge happens, just use original boundary data z
            bigz{j}=mergeZ{1}; bigz{j+1}=mergeZ{2};
        end
    end
end

if mcnt>imcnt % if fires have merged this step
    bigz = ROSAsmooth_v1_1(bigz,mcnt,J,resl,inswt,imswt); % smooth the fire lines
end
end

%% Appendix C1: Merge Function
function mergeZ = firemerge(z1,z2,bigz,j,J,shswt)
% = determine if fire lines z1 and z2 overlap and merge them with the
%   MATLAB "union" function if they do.
%
% Code:
x1 = real(z1); y1 = imag(z1); x1 = x1(~isnan(x1)); y1 = y1(~isnan(y1));
x2 = real(z2); y2 = imag(z2); x2 = x2(~isnan(x2)); y2 = y2(~isnan(y2)); % remove any NaN points from data.

poly1 = polyshape(x1,y1); poly2 = polyshape(x2,y2); % convert to polyshape objects.

if overlaps(poly1,poly2) % determine if the fire lines overlap
    polyun = union(poly1,poly2); [X,Y] = boundary(polyun); Z = X+1i*Y; % create union (merge) of the two lines and extract the boundary.
    if ~ispolycw(X,Y), Z = flip(Z); end % flip if anticlockwise - polygon must be oriented clockwise such that normal vector in correct direction.
    Pts = size(Z,1); smno=20; Zsm = firesmooth(Z,Pts,smno); % smooth the fire line.
    
    % Shift start/end of fire line away from fire junctions - areas of high curvature.
    Vertices = [real(Zsm.'); imag(Zsm.')].'; curvZ = LineCurvature2D(Vertices); % curvature of each point on the merged fire line.
    cmax = max(abs(curvZ)); pos = curvZ==cmax; % identifies the maximum curvature and its position in the boundary data.
    zmpk = Zsm(pos); % identifies the point with the highest curvature - the "peak point".
    zdis = abs(Zsm-zmpk.*ones(size(Zsm))); % find distance between "peak point" and all other boundary data.
    if J~=2 % if there are three fires, choose beginning point furthest from the other wildfire.
        if shswt == 124
            z3 = bigz{1}; % in the road example, fires 2 and 3 merge.
        else
            z3 = bigz{j+2}; % otherwise, assume fires 1 and 2 merge, leaving fire 3.
        end
        zdis2 = abs(Zsm-z3(1)*ones(size(Zsm)));
        zdis = zdis+zdis2; % adds the distance of each boundary point to the other wildfire.
    end
    zdmax = max(zdis); % finds distance of the point furthest from the peak point (plus the other wildfire, if applicable).
    
    while zdis(1)<zdmax % shifts Z coordinates until the start pt is the pt furthest from the peak point (and other wildfire, if applicable).
        Zsm = circshift(Zsm,-1); zdis = circshift(zdis,-1);
    end 
    mergeZ{1} = Zsm; % the final, new merged fire line.
else
    mergeZ{1} = z1; mergeZ{2} = z2; % output the input data if no overlap
end
end

%% Appendix C2: Fire Line Smooth and Interp Function
function znew = firesmooth(zold, Pts, smno)
% = smooths fire line and interpolates to desired resolution.
%   Uses "interparc" function [2].
%   Identical to Appendix B2 in ROSAsmooth.
%
% Code:
znew = zold;
for k = 1:smno % smooths a certain number of times.
    znew = smoothdata(znew,'gaussian', 4); znew = [znew; znew(1)]; % smooth and close the polygon.
    pt = interparc(Pts,real(znew), imag(znew),'spline'); % interpolates smoothed data.
    znew = pt(:,1)+1i*pt(:,2); % new boundary data.
end
end