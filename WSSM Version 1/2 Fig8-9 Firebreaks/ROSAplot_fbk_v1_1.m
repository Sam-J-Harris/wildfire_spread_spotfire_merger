function ROSAplot_fbk_v1_1(bigZ,bigJ,spc,shswt,imswt,fig)
% = plots the outputted fire spread. 
% Plots with a spacing 'spc' in figure number 'fig' (accounting for imswt).
% Code:
ssz = get(0,'ScreenSize'); set(gcf, 'Position',  [ssz(3)/4, ssz(4)/8, ssz(3)/2, ssz(3)/3])
Fsteps = size(bigJ,2); Ffig = fig + 2*imswt; % full number of steps (no spacing)

figure(Ffig)
set(gca,'XColor', 'none','YColor','none'), set(gca, 'color', 'none'); % transparent background and axes.
for k=1:spc:Fsteps % plot steps with spacing
    Jk = bigJ{k}; % number of fires J at step k
    for j=1:Jk
        z = bigZ{j,k}; x = real(z); y = imag(z);
        plot(x,y,'LineWidth',2), hold on,
    end
end
if shswt == 121 || shswt == 124 || shswt == 224
    xl = 4.5; Pts = 500; yln = linspace(-xl,xl,Pts);
    x1 = 1.5*ones(size(yln)); x2 = 1.7*ones(size(yln)); 
    
    yln3 = linspace(-xl,xl,200); x3 = 1.6*ones(size(yln3));
     
    lw = 2; plot(x1,yln,'k','LineWidth',lw); plot(x2,yln,'k','LineWidth',lw);
    patch([1.5 1.5 1.7 1.7], 1.1*[-xl xl xl -xl], 'k');
    plot(x3,yln3,'--w','LineWidth',2);

    set(gca,'XColor', 'none','YColor','none'), %set(gca, 'color', 'none');
elseif shswt == 122 % gap
    xl = 4; gp = 0.5; Pts = 500; 
    yln1 = linspace(-xl,-gp,ceil(Pts/2)); yln2 = linspace(gp,xl,ceil(Pts/2));
    x1 = 1.5*ones(size(yln1)); x2 = 1.7*ones(size(yln1));

    x3 = linspace(1.5,1.7,ceil(Pts/2)); x4 = linspace(1.5,1.7,ceil(Pts/2));
    y3 = -gp.*ones(size(x3)); y4 = gp.*ones(size(x4));
    
    lw = 2; plot(x1,yln1,'k','LineWidth',lw); plot(x2,yln1,'k','LineWidth',lw);
    plot(x1,yln2,'k','LineWidth',lw); 
    plot(x2,yln2,'k','LineWidth',lw);
    plot(x3,y3,'k','LineWidth',lw); 
    plot(x4,y4,'k','LineWidth',lw);
elseif shswt == 125 || shswt == 225 %lake
    xl = 5; lrad = 0.5; lz = lrad.*exp(1i*linspace(-pi,pi,500));
    lx = real(lz); ly = imag(lz);

    lw=2; %plot(lx,ly,'b','LineWidth',lw); 
    patch(lx,ly,'blue');

else
    xl = 3;
end
hold off, daspect([1 1 1]), axis([-xl xl -xl xl])
%title("Spread of multiple spotfires",'interpreter','latex','FontSize',18),

if shswt==1 || shswt==2 % additional plots for comparing with Hilton 2018 data
    scl = 1.7; axl = 5.0./scl; epsx = 0.0; epsy = 5.0./scl; % account for image scaling
    load('hilton_data.mat','I'); m=1; % loads background images

for k=1:spc:Fsteps
    figure(Ffig+m)
    set(gca,'XColor', 'none','YColor','none'), hold on, % transparent axes
    imagesc([-axl-epsx axl-epsx], [-2*axl+epsy epsy], I{m}), 
    set(gca,'YDir','reverse'), hold on, % scale image correctly and invert

    Jk = bigJ{k};
    for j=1:Jk
        z = bigZ{j,k}; x = real(z); y = imag(z);
        plot(x,-y,'y','LineWidth',3), 
    end
    daspect([1 1 1]); axis([-axl-epsx axl-epsx -(2*axl-epsy) epsy]);
    set(gca,'XColor', 'none','YColor','none') % transparent axes
    hold off

    m = m+1; if m>9, break, end % stop if run out of Hilton figures
end
hold off
end
end
