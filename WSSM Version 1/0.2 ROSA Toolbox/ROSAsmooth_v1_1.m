%% Appendix B: Fire Smooth Function
function bigz = ROSAsmooth_v1_1(bigz,mcnt,J,resl,inswt,imswt) % fire line smoothing
% ROSA smoothing function
%   smooths out instabilities in fire lines if inswt is turned on (=1).
%   Uses external functions "selfintersect" [1] and "interparc" [2]
%
% Refs:
%   [1] Canos, A.J., 2024. Fast and robust self-intersections. MATLAB 
%       Central File Exchange. Retrieved August 1, 2024. 
%       URL:(https://www.mathworks.com/matlabcentral/fileexchange/13351-fast-and-robust-self-intersections).
%
%   [2] D’Errico, J., 2024. interparc. MATLAB Central File Exchange. 
%       Retrieved August 1, 2024. 
%       URL: (https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc).
% 
% Code:
if inswt~=0 % interpolate and smooth to keep shapes at same resl, if inswt is switched on (=1).
    for j=1:J
        zv = bigz{j}; 
        if mcnt>0
            zv = fireintersect(zv,imswt); % checks fire lines don't intersect themselves (only perform after first merger)
        end
        if inswt==1
            polyw=polyshape(real(zv),imag(zv),'KeepCollinearPoints',true); Pts=round(resl.*perimeter(polyw)); 
        else
            Pts = size(zv,1); % keeps same number of pts for RK2 or RK4 timestepping
        end
        bigz{j} = firesmooth(zv,Pts,5); % smooth and interpolate fire line by desired no. of pts
    end
end
end

%% Appendix B1: Self-Intersect Check
function znew = fireintersect(z,imswt)
% = determine if a fire line intersects itself and removes any loops that do.
%   Uses "selfintersect" function [1].
%
% Code:
zt = z(1:end-1); xt = real(zt); yt=imag(zt);
[x0,y0,segments]=selfintersect(xt,yt); znew = z; % uses selfintersect function - see [1].

if size(segments,2)~=0 && size(segments,1)~=0 % ie there are overlapping segments
    SGstart = segments(:,1); SGend = segments(:,2); SGsize = size(SGstart,1); % determines number of overlapping segments
    for k=1:SGsize
        zt(SGstart(k))=x0(k)+1i*y0(k); % retain the first point of the segment
        zt(SGstart(k)+1:SGend(k))=0;
        zt(SGstart(k):SGend(k))=0; % turn all other points in the segment to zero
    end
    zt(zt==0)=[]; Pts = size(zt,1); smno = 1;
    znew = firesmooth(zt,Pts,smno); % remove all 0 points from fire line data and smooths curve
    ROSAsimage_v1_1(z,zt,x0,y0,imswt) % shows image of overlapping segments,  if imswt is turned on.
end
end

%% Appendix B2: Fire Line Smooth and Interp Function
function znew = firesmooth(zold, Pts, smno)
% = smooths fire line and interpolates to desired resolution.
%   Uses "interparc" function [2].
%   Identical to Appendix C2 in ROSAmerger.
%
% Code:
znew = zold;
for k = 1:smno % smooths a certain number of times.
    znew = smoothdata(znew,'gaussian', 4); znew = [znew; znew(1)]; % smooth and close the polygon.
    pt = interparc(Pts,real(znew), imag(znew),'spline'); % interpolates smoothed data.
    znew = pt(:,1)+1i*pt(:,2); % new boundary data.
end
end