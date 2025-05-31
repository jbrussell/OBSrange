%% plot the misfit histogram
f2 = figure(2); clf;
[ Ncount_Erms, cent_Erms ] = plot_hist(gca,E_rms*1000,Nbins);
set(gca,'fontsize',16,'linewidth',2); box on;
title('\textbf{Misfit}','fontsize',18,'interpreter','latex');
xlabel('RMS (ms)','fontsize',18,'interpreter','latex');

% Plot value for full dataset (bootstrap index = 1);
xline(gca,E_rms(1)*1000,'-','color',[147, 233, 190]/255,'alpha',1,'linewidth',3);
