function ROSAsimage_v1_1(z,zt,x0,y0,imswt) 
% = plot of overlapping segments of a fire line curve
if imswt==1 % only if imswt turned on (=1).
    figure(2)
    scsz = get(0,'ScreenSize'); 
    set(gcf, 'Position',  [scsz(3)/1.8, scsz(4)/8, scsz(3)/2.5, scsz(3)/4]) % position based on user's screensize
    plot(real(z),imag(z),'b',x0,y0,'.r'), % z = full curve, x0,y0 = point of intersection(s) 
    hold on, plot(real(zt),imag(zt),'g'), hold off, % zt = non-overlapping final curve
    axis ('equal'); grid
end
end