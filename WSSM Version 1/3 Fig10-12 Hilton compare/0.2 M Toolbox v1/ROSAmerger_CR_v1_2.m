%% Appendix B: Fire Merger function
function [bigz, bigc, merdata, J] = ROSAmerger_CR_v1_2(bigz,bigc,merdata,J,resl,intswit,imswit)
    mswit = 0; % variable to assess if fires have merged this step
    if J~=1
        mergeZ = firemerge(bigz{1}, bigz{2});
        if size(mergeZ,2)==1
            premerge = bigz; % original z data (two fires)
            bigz{1} = mergeZ{1}; J=1; mswit=1; % merged fire

            zcf = mergeZ{1}; xcf = real(zcf); ycf = imag(zcf); % find centre of new polygon
            polycf = polyshape(xcf,ycf); [xc, yc] = centroid(polycf); 
            ctr = xc+1i*yc; ctr = -1;
            bigc = {}; bigc{1} = ctr; % new centre
            
            postmerge = bigz; % merged z data (one fire)
            merdata = {premerge, postmerge};
        else
            bigz{1}=mergeZ{1}; bigz{2}=mergeZ{2};
        end
    else
        bigz{1} = fireintersect(bigz{1},imswit);
    end

    % Fire Interpolation
    if intswit==1 && mswit==0 %interpolate to keep shapes at same resl#
        smno = 2;
        for j=1:J
            zv = bigz{j}; xv=real(zv); yv=imag(zv); 
            polyw=polyshape(xv,yv,'KeepCollinearPoints',true); Pts=round(resl.*perimeter(polyw)); 
            bigz{j} = firesmooth(zv,Pts,smno);
        end
    end
end

%% Appendix B1: Merge Function
function mergeZ = firemerge(z1,z2)
    x1 = real(z1); y1 = imag(z1); x1 = x1(~isnan(x1)); y1 = y1(~isnan(y1));
    x2 = real(z2); y2 = imag(z2); x2 = x2(~isnan(x2)); y2 = y2(~isnan(y2));

    poly1 = polyshape(x1,y1); poly2 = polyshape(x2,y2);

    if overlaps(poly1,poly2)
        polyun = union(poly1,poly2); [X,Y] = boundary(polyun); Z = X+1i*Y; 
        if ~ispolycw(X,Y), Z = flip(Z); end % flip if anticlockwise
        Pts = size(Z,1); smno=20; 
        Zsm = firesmooth(Z,Pts,smno); Zshf = circshift(Zsm,round(Pts/2)); 
        mergeZ{1} = Zshf;
    else
        mergeZ{1} = z1; mergeZ{2} = z2;
    end
end

%% Appendix B2: Self-Intersect Check
function znew = fireintersect(z,imswit)
    zt = z(1:end-1); xt = real(zt); yt=imag(zt);
    [x0,y0,segments]=selfintersect(xt,yt); znew = z; 
    
    if size(segments,2)~=0
        SGstart = segments(:,1); SGend = segments(:,2); SGsize = size(SGstart,1);
        for k=1:SGsize
            zt(SGstart(k))=x0(k)+1i*y0(k);
            zt(SGstart(k)+1:SGend(k))=0;
            zt(SGstart(k):SGend(k))=0;
        end
        zt(zt==0)=[]; Pts = size(zt,1); smno = 1;
        znew = firesmooth(zt,Pts,smno);
    
        if imswit==1
            figure(2)
            scsz = get(0,'ScreenSize'); 
            set(gcf, 'Position',  [scsz(3)/1.8, scsz(4)/8, scsz(3)/2.5, scsz(3)/4])
            plot(real(z),imag(z),'b',x0,y0,'.r'), hold on, 
            plot(real(zt),imag(zt),'g'), hold off, axis ('equal'); grid
        end
    end

end

%% Appendix B3: Order Points Around Polygon
function znew = fireorder(zold)
    zuni = unique(zold); oldlst = zuni(~isnan(zuni));
    M = size(oldlst,1); 
    newlst = [oldlst(1)]; oldlst = oldlst(oldlst~=oldlst(1));
    for j=1:M
        zcurr = newlst(end); Mnew = size(oldlst,1); DisMat = ones(Mnew,1);
        for k=1:Mnew
            DisMat(k) = norm(zcurr-oldlst(k));
        end
        minDis = min(DisMat);
        znext = oldlst(DisMat==minDis); 
        if minDis<0.1
            newlst = [newlst znext]; 
        end    
        oldlst = oldlst(oldlst~=znext);
    end
    Mfin = size(newlst,1);
    znew = circshift(newlst.', round(Mfin/4));
    if ~ispolycw(real(znew),imag(znew)), znew = flip(znew); end
end

%% Appendix B4: Fire Line Smooth and Interp Function
function znew = firesmooth(zold, Pts, smno)
    znew = zold;
    for k = 1:smno
        znew = smoothdata(znew,'gaussian',4); znew = [znew; znew(1)];
        znew = fireinterp(znew,Pts,'spline'); 
    end
end