function ROSSOplotRE_v1_5(tvec,bigRE,fignum)
scrsiz = get(0,'ScreenSize'); set(gcf, 'Position',  [scrsiz(3)/4, scrsiz(4)/8, scrsiz(3)/2.5, scrsiz(3)/4])

figure(fignum)
plot(tvec,bigRE{1},tvec,bigRE{2},tvec,bigRE{3},'LineWidth',1.25), 
legend({'N=32','N=64','N=128'},'location','northwest','FontSize',25), 
axis square; set(gca,'FontSize',18)
end