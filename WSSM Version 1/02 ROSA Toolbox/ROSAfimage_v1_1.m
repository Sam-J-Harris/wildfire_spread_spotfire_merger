function ROSAfimage_v1_1(bigZ,bigJ,m,imswt) 
% = plot of fire line evolution during the ROSA algorithm
if imswt==1 % only if imswt turned on (=1).
    figure(1) % add to figure 1 on each time step
    scsz = get(0,'ScreenSize'); set(gcf, 'Position',  [scsz(3)/20, scsz(4)/8, scsz(3)/2, scsz(3)/3])
    axl=4; axes('position',[.02 .3 .5 .5]);
    for k=1:m+1
        for j=1:bigJ{k} % plot all j fire lines on each time step from 1 to m+1
            plot(real(bigZ{j,k}),imag(bigZ{j,k})), hold on, 
        end
        if m==1 && k==1
            axis([-axl axl -axl axl]), daspect([1 1 1]), pause(0.05)
        end
    end
    hold off, axis([-axl axl -axl axl]), axis square, %daspect([1 1 1]),
    title("Wildfires",'interpreter','latex'),
end
end