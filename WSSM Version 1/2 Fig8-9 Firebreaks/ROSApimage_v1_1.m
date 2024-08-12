function ROSApimage_v1_1(Allpols,bigz,J,phi,imswt) % phi contour image
    inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % ftn to determine if point is inside polygon
    if imswt==1
        axl = 7; LW = 'linewidth'; MS = 'markersize'; ms = 6; PO = 'position';
        figure(1)
        axes(PO,[.5 .3 .5 .5]);
        plot(Allpols,'.r',MS,ms), hold on
        x = linspace(-axl,axl,500); [xx,yy] = meshgrid(x,x); zz = xx+1i*yy; 
        for j=1:J
            if j==1
                cond = inpolygonc(zz,bigz{j});
            else
                cond = cond|inpolygonc(zz,bigz{j});
            end
        end
        uu = phi(zz); uu(cond) = NaN;
        for j=1:J % plotting boundaries
            zb = bigz{j}; xb = real(zb); yb = imag(zb); plot(xb,yb,'k',LW,.9), 
        end
        contourf(x,x,uu,30,LW,1), hold off, 
        axis([-axl axl -axl axl]), axis square, title('Laplace solution');
    end
end