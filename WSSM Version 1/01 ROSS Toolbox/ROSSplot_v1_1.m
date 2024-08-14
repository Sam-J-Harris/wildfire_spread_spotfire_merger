function ROSSplot_v1_1(bigZ,spc,fignum)
% Plots the wildfire evolution
%
% Inputs:
%   bigZ = {Z1,Z2,...Zn}: list of n experiments of the wildfire boundary data.
%   spc = determines how many isochrones are plotted - set spc = 1 to plot all.
%   fignum = which figure the plots are presented in.
%
% Code:
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4]) % outputs the plots in specific position, dependent on user's screensize

Zno = size(bigZ,2); % number of experiments to be plotted, these appear in one long row
for j=1:Zno
    Z = bigZ{j}; % jth experiment of a single fire evolution
    ftstep = size(Z,2); %final # of steps - note that ftstep =< itstep
    
    figure(fignum)
    subplot(1, Zno, j); % puts fire evolution in the jth subfigure
    set(gca,'XColor', 'none','YColor','none')
    hold on %uncomment/comment "hold" if you want/do not want to see each level set
    for k=1:spc:ftstep
        z = Z{k}; x = real(z); y=imag(z);
    
        figure(fignum)
        subplot(1, Zno, j);
        plot(x,y,'LineWidth',1.25), axis square; daspect([1 1 1]);
    end
end
end