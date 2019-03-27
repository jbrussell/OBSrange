%% plot the misfit histogram
f2 = figure(2); clf;
[ Ncount_Erms, cent_Erms ] = plot_hist(gca,E_rms*1000,Nbins);
set(gca,'fontsize',16,'linewidth',2); box on;
title('\textbf{Misfit}','fontsize',18,'interpreter','latex');
xlabel('RMS (ms)','fontsize',18,'interpreter','latex');
