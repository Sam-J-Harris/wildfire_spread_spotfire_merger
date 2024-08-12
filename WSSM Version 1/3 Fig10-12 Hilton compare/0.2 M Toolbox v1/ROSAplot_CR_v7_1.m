function ROSAplot_CR_v7_1(bigZ,bigJ,spc,shswit,fignum)
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])
fsteps = size(bigJ,2);

figure(fignum)
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
for k=1:spc:fsteps
    Jk = bigJ{k};
    for j=1:Jk
        z = bigZ{j,k}; x = real(z); y = imag(z);
        plot(x,y,'LineWidth',1.25), hold on,
    end
end
hold off,  axis square, axis([-3 3 -3 3]);
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])
%title("Merging of two circular wildfires",'interpreter','latex','FontSize',18),
%set(gca,'XColor', 'none','YColor','none'), hold on, %set(gca, 'color', 'none');
% ax1 = gca;                   % gca = get current axis
% ax1.YAxis.Visible = 'off';   % remove y-axis
% ax1.XAxis.Visible = 'off';   % remove x-axis
set(gca, 'YTick', []); set(gca, 'XTick', []);

if shswit~=0
    scl = 1.7;
    I = {}; axl = 5.0./scl; epsx = 0.0; epsy = 5.0./scl; 
    I{1} = imread("hilton2018_fig8_05_csi.png");
    I{2} = imread("hilton2018_fig8_10_csi.png");
    I{3} = imread("hilton2018_fig8_15_csi.png"); 
    I{4} = imread("hilton2018_fig8_20_csi.png"); 
    I{5} = imread("hilton2018_fig8_25_csi.png");
    I{6} = imread("hilton2018_fig8_30_csi.png");
    I{7} = imread("hilton2018_fig8_35_csi.png");
    I{8} = imread("hilton2018_fig8_40_csi.png");
    I{9} = imread("hilton2018_fig8_45_csi.png");

m=1;
for k=1:spc:fsteps

    figure(fignum+m)
    set(gca,'XColor', 'none','YColor','none'), hold on, %set(gca, 'color', 'none');
    imagesc([-axl-epsx axl-epsx], [-2*axl+epsy epsy], I{m}), 
    set(gca,'YDir','reverse'), hold on,

    Jk = bigJ{k};
    for j=1:Jk
        z = bigZ{j,k}; x = real(z); y = imag(z);
        plot(x,-y,'y','LineWidth',3), 
    end
    daspect([1 1 1]); axis([-axl-epsx axl-epsx -(2*axl-epsy) epsy]);
    set(gca,'XColor', 'none','YColor','none') %set(gca, 'color', 'none');
    scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])
    hold off

    m = m+1;
end
hold off
end
end
