function ROSSOplotPY_v1_5(bigZ,spc,fignum)
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])

Zno = size(bigZ,2);
for j=1:Zno
    Z = bigZ{j};
    ftstep = size(Z,2); %final # of steps - note that ftstep =< itstep
    
    figure(fignum)
    subplot(1, Zno, j);
    set(gca,'XColor', 'none','YColor','none') %set(gca, 'color', 'none');
    hold on %uncomment/comment "hold" if you want/do not want to see each level set
    for k=1:spc:ftstep
    z = Z{k}; x = real(z); y=imag(z);

    figure(fignum)
    subplot(1, Zno, j);
    plot(x,y,'LineWidth',1.25), axis([-3 3 -1 6]); daspect([1 1 1]);
    end
end
end

