clear all
load('ap_removed.mat');

if length(frex) == 391
    frex = frex(21:end);
end
mean_spect = ap_removed;

delta = [5 6]; delta_idx = frex >= delta(1) & frex < delta(2);
theta = [6 8]; theta_idx = frex >= theta(1) & frex < theta(2);
alpha = [8 14]; alpha_idx = frex >= alpha(1) & frex < alpha(2);
beta = [14 30]; beta_idx = frex >= beta(1) & frex < beta(2);

delta_power_ap = trimmean(mean_spect(:,:,delta_idx),0,'round',3);
theta_power_ap = trimmean(mean_spect(:,:,theta_idx),0,'round',3);
alpha_power_ap = trimmean(mean_spect(:,:,alpha_idx),0,'round',3);
beta_power_ap = trimmean(mean_spect(:,:,beta_idx),0,'round',3);

delta_mean = mean(delta_power_ap, 1);
theta_mean = mean(theta_power_ap, 1);
alpha_mean = mean(alpha_power_ap, 1);
beta_mean = mean(beta_power_ap, 1);

delta_sem = withinstde(delta_power_ap);
theta_sem = withinstde(theta_power_ap);
alpha_sem = withinstde(alpha_power_ap);
beta_sem = withinstde(beta_power_ap);

SPR_ap = (alpha_power_ap + beta_power_ap) ./ (delta_power_ap + theta_power_ap);
SPR_mean = mean(SPR_ap,1);
SPR_sem = withinstde(SPR_ap);

[~,p_12,~,stat_12] = ttest(SPR_ap(:,2), SPR_ap(:,1));
[~,p_13,~,stat_13] = ttest(SPR_ap(:,3), SPR_ap(:,1));
[~,p_23,~,stat_23] = ttest(SPR_ap(:,3), SPR_ap(:,2));


figure('Position',[76, 76, 870, 540]); 
subplot(2,4,1); hold on;
plot(1:3, delta_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, delta_mean, delta_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)');
ylim([-.03 .6]);
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'', '', ''});
title('Delta')

subplot(2,4,2); hold on;
plot(1:3, theta_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, theta_mean, theta_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)')
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'', '', ''});
title('Theta')

subplot(2,4,5); hold on;
plot(1:3, alpha_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, alpha_mean, alpha_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)');
ylim([-.12 .7]);
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'Baseline', 'Follow up', 'End of study'});
title('Alpha')

subplot(2,4,6); hold on;
plot(1:3, beta_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, beta_mean, beta_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)');
ylim([-.1 .3]);
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'Baseline', 'Follow up', 'End of study'});
title('Beta')

subplot(1,2,2);hold on;
plot(1:3, SPR_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, SPR_mean, SPR_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
title('Spectral Power Ratio');
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'Baseline', 'Follow up', 'End of study'});


%%
clear all
load('mean_spect.mat');

delta = [3 4]; delta_idx = frex >= delta(1) & frex < delta(2);
theta = [4 8]; theta_idx = frex >= theta(1) & frex < theta(2);
alpha = [8 14]; alpha_idx = frex >= alpha(1) & frex < alpha(2);
beta = [14 30]; beta_idx = frex >= beta(1) & frex < beta(2);

delta_power_ap = trimmean(mean_spect(:,:,delta_idx),0,'round',3);
theta_power_ap = trimmean(mean_spect(:,:,theta_idx),0,'round',3);
alpha_power_ap = trimmean(mean_spect(:,:,alpha_idx),0,'round',3);
beta_power_ap = trimmean(mean_spect(:,:,beta_idx),0,'round',3);

delta_mean = mean(delta_power_ap, 1);
theta_mean = mean(theta_power_ap, 1);
alpha_mean = mean(alpha_power_ap, 1);
beta_mean = mean(beta_power_ap, 1);

delta_sem = withinstde(delta_power_ap);
theta_sem = withinstde(theta_power_ap);
alpha_sem = withinstde(alpha_power_ap);
beta_sem = withinstde(beta_power_ap);

SPR_ap = (alpha_power_ap + beta_power_ap) ./ (delta_power_ap + theta_power_ap);
SPR_mean = mean(SPR_ap,1);
SPR_sem = withinstde(SPR_ap);

figure('Position',[76, 76, 870, 540]); 
subplot(2,4,1); hold on;
plot(1:3, delta_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, delta_mean, delta_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)');
% ylim([-.03 .6]);
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'', '', ''});
title('Delta')

subplot(2,4,2); hold on;
plot(1:3, theta_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, theta_mean, theta_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)')
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'', '', ''});
title('Theta')

subplot(2,4,5); hold on;
plot(1:3, alpha_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, alpha_mean, alpha_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)');
% ylim([-.12 .7]);
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'Baseline', 'Follow up', 'End of study'});
title('Alpha')

subplot(2,4,6); hold on;
plot(1:3, beta_power_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, beta_mean, beta_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
ylabel('Power (\muV^2)');
ylim([0 .25]);
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'Baseline', 'Follow up', 'End of study'});
title('Beta')

subplot(1,2,2);hold on;
plot(1:3, SPR_ap, 'linewidth', 1, 'Color', [.7 .7 .7]);
errorbar(1:3, SPR_mean, SPR_sem, 'color', [0 0 0], 'linewidth', 1.5);
set(gca,'linewidth', 1, 'fontsize', 11, 'fontweight','bold', 'box', 'off');
title('Spectral Power Ratio');
xlim([.8 3.2]);
xticks([1 2 3]);
% xticklabels({'BL', 'FU', 'EOS'});
xticklabels({'Baseline', 'Follow up', 'End of study'});


