function ROSApimage_v1_1(Allpols,bigz,J,phi,imswt) % phi contour image
% = contour plot of exterior phi calculated from AAA-LS algorithm.
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % ftn to determine if point is inside polygon
if imswt==1 % only if imswt turned on (=1).
    axl = 7; LW = 'linewidth'; MS = 'markersize'; ms = 6; PO = 'position';
    figure(1)
    axes(PO,[.5 .3 .5 .5]);
    plot(Allpols,'.r',MS,ms), hold on % plot of interior poles found.
    x = linspace(-axl,axl,500); [xx,yy] = meshgrid(x,x); zz = xx+1i*yy; 
    for j=1:J % only plot contours outside of ALL fire line curves
        if j==1
            cond = inpolygonc(zz,bigz{j});
        else
            cond = cond|inpolygonc(zz,bigz{j}); 
        end
    end
    uu = phi(zz); uu(cond) = NaN;
    for j=1:J % plotting fire lines
        zb = bigz{j}; xb = real(zb); yb = imag(zb); plot(xb,yb,'k',LW,.9), 
    end
    contourf(x,x,uu,30,LW,1), hold off, % filled contour plot
    axis([-axl axl -axl axl]), axis square, title('Laplace solution');
end
end