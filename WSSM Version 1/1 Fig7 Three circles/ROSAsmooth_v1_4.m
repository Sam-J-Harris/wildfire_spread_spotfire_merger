%% Appendix B: Fire Smooth Function
function bigz = ROSAsmooth_v1_4(bigz,mcnt,J,resl,inswt,imswt) % fire line smoothing
if inswt~=0 %interpolate and smooth to keep shapes at same resl
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
% Uses selfintersect function from [REF online].
% Code:
zt = z(1:end-1); xt = real(zt); yt=imag(zt);
[x0,y0,segments]=selfintersect(xt,yt); znew = z; % uses selfintersect function - see REF.

if size(segments,2)~=0 && size(segments,1)~=0 % ie there are overlapping segments
    SGstart = segments(:,1); SGend = segments(:,2); SGsize = size(SGstart,1); % determines number of overlapping segments
    for k=1:SGsize
        zt(SGstart(k))=x0(k)+1i*y0(k); % retain the first point of the segment
        zt(SGstart(k)+1:SGend(k))=0;
        zt(SGstart(k):SGend(k))=0; % turn all other points in the segment to zero
    end
    zt(zt==0)=[]; Pts = size(zt,1); smno = 1;
    znew = firesmooth(zt,Pts,smno); % remove all 0 points from fire line data and smooth curve
    ROSAsimage_v1_1(z,zt,x0,y0,imswt) % shows image of overlapping segments, if selected.
end
end

%% Appendix B2: Fire Line Smooth and Interp Function
function znew = firesmooth(zold, Pts, smno)
% = smooths fire line and interpolates to desired numbers of pts
% Code:
znew = zold;
for k = 1:smno % smooths a certain number of times
    znew = smoothdata(znew,'gaussian', 4); znew = [znew; znew(1)]; % smooth and close
    znew = fireinterp(znew,Pts,'spline'); % interpolate smoothed data
end
end