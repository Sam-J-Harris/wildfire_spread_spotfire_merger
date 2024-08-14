function ROSAplot_v1_1(bigZ,bigJ,spc,shswt,imswt,fig)
% = plots the final, outputted evolution of the fire lines with a spacing 
%   'spc' in figure number 'fig' (accounting for if imswt was used 
%   during ROSA).
%
% Code:
ssz = get(0,'ScreenSize'); set(gcf, 'Position',  [ssz(3)/4, ssz(4)/8, ssz(3)/2, ssz(3)/3]) % gets user's screensize.
Fsteps = size(bigJ,2); Ffig = fig + 2*imswt; % full number of steps (no spacing); figure number accounting for imswt.

figure(Ffig)
set(gca,'XColor', 'none','YColor','none'), set(gca, 'color', 'none'); % transparent background and axes.
for k=1:spc:Fsteps % plot steps with spacing.
    Jk = bigJ{k}; % number of fires at step k
    for j=1:Jk
        z = bigZ{j,k}; x = real(z); y = imag(z); % plots all fire lines at step k.
        plot(x,y,'LineWidth',1.25), hold on,
    end
end
hold off, daspect([1 1 1]),
if shswt==07, axis([-3 3 -3 3]); end % fix axes for three circles example.
title("Spread of multiple spotfires",'interpreter','latex','FontSize',18),

%% Hilton et al. 2018 additional figures
%   The following can only be used if the user has access to the experimental
%   images from Hilton et al. (2018)/ Sullivan et al. (2019). 
% 
%   If so, update the "load background image" line to refelct which folder 
%   they are in, then uncomment the following code.
%   
%   If not, the above figure will be just fine.
%
% Code:
% if shswt==1 || shswt==2 % additional plots for comparing with Hilton 2018 data
%     scl = 1.7; axl = 5.0./scl; epsx = 0.0; epsy = 5.0./scl; % account for image scaling
%     load('hilton_data.mat','I'); m=1; % loads background images
% 
% for k=1:spc:Fsteps
%     figure(Ffig+m)
%     set(gca,'XColor', 'none','YColor','none'), hold on, % transparent axes
%     imagesc([-axl-epsx axl-epsx], [-2*axl+epsy epsy], I{m}), 
%     set(gca,'YDir','reverse'), hold on, % scale image correctly and invert
% 
%     Jk = bigJ{k};
%     for j=1:Jk
%         z = bigZ{j,k}; x = real(z); y = imag(z);
%         plot(x,-y,'y','LineWidth',3), 
%     end
%     daspect([1 1 1]); axis([-axl-epsx axl-epsx -(2*axl-epsy) epsy]);
%     set(gca,'XColor', 'none','YColor','none') % transparent axes
%     hold off
% 
%     m = m+1; if m>9, break, end % stop if run out of Hilton figures
% end
% hold off
% end
end
