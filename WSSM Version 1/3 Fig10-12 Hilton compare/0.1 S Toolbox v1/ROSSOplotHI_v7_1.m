function ROSSOplotHI_v7_1(Z,input,fignum)
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])
ftstep = size(Z,2); %final # of steps - note that ftstep =< itstep
I={[],[],[],[]}; axl = 1.95; epsx = 0.025; epsy = 0.75;

if input==185
    load('scl_hilton2018_fig5b_v4.mat','scl');
    I = {}; axl = 5.0./scl; epsx = 0.0; epsy = 5.0./scl;
    I{1} = imread("hilton2018_fig5b.png");
    I{2} = imread("hilton2018_fig5c.png"); 
    I{3} = imread("hilton2018_fig5d.png"); 
    I{4} = imread("hilton2018_fig5e.png"); 
elseif input==187
    load('scl_hilton2018_fig7e_csi_v4.mat','scl');
    I = {}; axl = 5.0./scl; epsx = 0.0; epsy = 5.0./scl; 
    I{1} = imread("hilton2018_fig7e_csi.png");
    I{2} = imread("hilton2018_fig7f_csi.png"); 
    I{3} = imread("hilton2018_fig7g_csi.png"); 
    I{4} = imread("hilton2018_fig7h_csi.png"); 
end

figure(fignum)
set(gca,'XColor', 'none','YColor','none') %set(gca, 'color', 'none');
hold on %uncomment/comment "hold" if you want/do not want to see each level set
for k=1:ftstep
    z = Z{k}; x = real(z); y=imag(z);

    figure(fignum)
    subplot(2,2,k)
    imagesc([-axl-epsx axl-epsx], [-2*axl+epsy epsy], I{k}), 
    set(gca,'YDir','reverse'), hold on,
    plot(x,-y,'y','LineWidth',3), 
    daspect([1 1 1]); axis([-axl-epsx axl-epsx -(2*axl-epsy) epsy]);
    set(gca,'XColor', 'none','YColor','none') %set(gca, 'color', 'none');
end
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])
hold off
end

