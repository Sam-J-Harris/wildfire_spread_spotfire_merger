function ROSSplotRE_v1_1(tvec,bigRE,fignum)
% Plots the relative error of the rate of change of area law over time
% 
% Inputs:
%   tvec = time vector.
%   bigRE = {RE1, RE2,...,REn} = relative error over time of n wildfire experiments.
%   fignum = which figure the plots are presented in.
%
% Code:
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4]) % outputs the plots in specific position, dependent on user's screensize

figure(fignum)
plot(tvec,bigRE{1},tvec,bigRE{2},tvec,bigRE{3},'LineWidth',1.25), % plots of relative error vs time
legend({'N=32','N=64','N=128'},'location','northwest','FontSize',25), % edit legend as required
axis square; set(gca,'FontSize',18)
end