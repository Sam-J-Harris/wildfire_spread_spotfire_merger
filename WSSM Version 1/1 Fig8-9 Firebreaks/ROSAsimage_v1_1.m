function ROSAsimage_v1_1(z,zt,x0,y0,imswt) % segment removal image
if imswt==1
    figure(2)
    scsz = get(0,'ScreenSize'); 
    set(gcf, 'Position',  [scsz(3)/1.8, scsz(4)/8, scsz(3)/2.5, scsz(3)/4])
    plot(real(z),imag(z),'b',x0,y0,'.r'), 
    hold on, plot(real(zt),imag(zt),'g'), hold off, 
    axis ('equal'); grid
end
end