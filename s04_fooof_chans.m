% 
folders = dir('V:\EEG\*');
folders = folders(3:end);

num_sessions = 3;
sess_ids = {'BL','FU','EOS'};
num_ps = 12;

frex = 3:.1:40;

load('mean_spect.mat');
% 
py.importlib.import_module('os');
py.importlib.import_module('numpy');
py.importlib.import_module('fooof');

settings = struct(); 
settings.peak_width_limits = [.5 9];
settings.max_n_peaks = 6;
settings.min_peak_height = 0;
settings.peak_threshold = 2.0;
settings.aperiodic_mode = 'fixed';
settings.verbose = true;
f_range = [3 40];

% ap_removed = nan(size(mean_spect));
slope = nan(num_ps,num_sessions,64);
offset = slope;

for p = 1:num_ps

    for s = 1:num_sessions
        session_id = sess_ids{s};
        location = [folders(p).folder '\' folders(p).name '\' folders(p).name '_' session_id '\' folders(p).name '_' session_id];
        load([location '_spect.mat'],'spect');

        for c = 1:64

            fooof_results(p,s,c) = fooof(frex, squeeze(spect(c,:)), f_range, settings, true); %#ok<SAGROW> 
            ap_removed(p,s,c,:) = fooof_results(p,s,c).power_spectrum - fooof_results(p,s,c).ap_fit; %#ok<SAGROW> 
            total_spect(p,s,c,:) = fooof_results(p,s,c).power_spectrum; %#ok<SAGROW> 

            slope(p,s,c) = fooof_results(p,s,c).aperiodic_params(2);
            offset(p,s,c) = fooof_results(p,s,c).aperiodic_params(1);
        end
    end
end

% %save fooof params and residual spectra separately
save('fooof_results_chan.mat','fooof_results', 'slope', 'offset');
save('ap_removed_chan.mat','ap_removed','total_spect','frex');

return
%%
load('fooof_results_chan.mat');
load('ap_removed_chan.mat');

load('chans.mat');
figure; 
subplot(231);
topoplot(squeeze(mean(slope(:,1,:),1)),chans(1:64),'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
caxis([.8 1.8]);
colorbar;
subplot(232);
topoplot(squeeze(mean(slope(:,2,:),1)),chans(1:64),'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
caxis([.8 1.8]);
colorbar;
subplot(233);
topoplot(squeeze(mean(slope(:,3,:),1)),chans(1:64),'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
caxis([.8 1.8]);
colorbar;

subplot(234);
topoplot(squeeze(mean(slope(:,2,:) - slope(:,1,:),1)),chans(1:64),'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
caxis([-.4 .4]);
colorbar;
title('1 v 2');
subplot(235);
topoplot(squeeze(mean(slope(:,3,:) - slope(:,2,:),1)),chans(1:64),'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
caxis([-.4 .4]);
colorbar;
title('2 v 3');

%%
[~,p12] = ttest(squeeze(slope(:,1,:)),squeeze(slope(:,2,:)));
[~,p13] = ttest(squeeze(slope(:,1,:)),squeeze(slope(:,3,:)));
[~,p23] = ttest(squeeze(slope(:,2,:)),squeeze(slope(:,3,:)));

nei = spatial_neighbors(chans(1:64),56);

[p_pos12, chans_pos12, p_neg12, chans_neg12] = ClustPerm1D(p12', nei, mean(squeeze(slope(:,1,:)) - squeeze(slope(:,2,:)),1)');
[p_pos13, chans_pos13, p_neg13, chans_neg13] = ClustPerm1D(p13', nei, mean(squeeze(slope(:,1,:)) - squeeze(slope(:,3,:)),1)');
[p_pos23, chans_pos23, p_neg23, chans_neg23] = ClustPerm1D(p23', nei, mean(squeeze(slope(:,2,:)) - squeeze(slope(:,3,:)),1)');

p12_corr = [p_pos12 p_neg12];
p13_corr = [p_pos13 p_neg13];
p23_corr = [p_pos23 p_neg23];

disp([num2str(min(p12)) ' ' num2str(min(p12_corr))]);
disp([num2str(min(p13)) ' ' num2str(min(p13_corr))]);
disp([num2str(min(p23)) ' ' num2str(min(p23_corr))]);

save('fooof_chans_statistics.mat');

%%

load('chans.mat');
figure; 
subplot(231);
topoplot(squeeze(mean(slope(:,1,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off', 'whitebk','on');
caxis([.8 1.8]);
h = colorbar;
title('Session 1');
ylabel(h,'Slope');
set(h,'fontsize',12,'fontweight','bold','linewidth',1)
set(gca,'fontsize',12,'fontweight','bold','linewidth',1)
subplot(232);
topoplot(squeeze(mean(slope(:,2,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off', 'whitebk','on');
caxis([.8 1.8]);
h = colorbar;
title('Session 2');
ylabel(h,'Slope');
set(h,'fontsize',12,'fontweight','bold','linewidth',1)
set(gca,'fontsize',12,'fontweight','bold','linewidth',1)
subplot(233);
topoplot(squeeze(mean(slope(:,3,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([.8 1.8]);
h = colorbar;
title('Session 3');
ylabel(h,'Slope');
set(h,'fontsize',12,'fontweight','bold','linewidth',1)
set(gca,'fontsize',12,'fontweight','bold','linewidth',1)
% figure;
subplot(234);
topoplot(squeeze(mean(slope(:,2,:) - slope(:,1,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off', ...
    'emarker2',{find(chans_neg12), '.', 'w', 15},'whitebk','on');
caxis([-.4 .4]);
h = colorbar;
title('Session 2 vs 1');
ylabel(h,'Slope difference');
set(h,'fontsize',12,'fontweight','bold','linewidth',1)
set(gca,'fontsize',12,'fontweight','bold','linewidth',1)
% figure;
subplot(235);
topoplot(squeeze(mean(slope(:,3,:) - slope(:,2,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([-.4 .4]);
h = colorbar;
title('Session 3 vs 2');
ylabel(h,'Slope difference');
set(h,'fontsize',12,'fontweight','bold','linewidth',1)
set(gca,'fontsize',12,'fontweight','bold','linewidth',1)
% figure;
subplot(236);
topoplot(squeeze(mean(slope(:,3,:) - slope(:,1,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off', ...
    'emarker2',{find(chans_neg13), '.', 'w', 15},'whitebk','on');
caxis([-.4 .4]);
h = colorbar;
title('Session 3 vs 1');
ylabel(h,'Slope difference');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);

%%

for p = 1:12
    for s = 1:3
        for c = 1:64
            total_spect(p,s,c,:) = fooof_results(p,s,c).power_spectrum; %#ok<SAGROW> 
        end
    end
end

%%
figure; 
loglog(frex, squeeze(mean(total_spect(:,:,find(chans_neg13),:),[1 3]))','linewidth',1.5);
legend('BL','FU','EOS'); legend boxoff;
xlabel('log10 Frequency (Hz)');
ylabel('log10 Power');
set(gca, 'linewidth',1,'fontsize',12,'fontweight','bold','box','off');
xlim([2.9 45]);
ylim([-2.1 -.15]);


%%
col = .75;

ms13 = mean(slope(:,:,find(chans_neg13)),3);
ms12 = mean(slope(:,:,find(chans_neg12)),3);

figure;
subplot(121);
hold on;
% plot(1:3, log10(SPR_ap), 'linewidth',1,'color',[col col col]);
% plot(1:3, mean(log10(SPR_ap),1),'linewidth',2,'color','k');
plot(1:3, ms12, 'linewidth',1,'color',[col col col]);
errorbar(1:3, mean(ms12,1),withinstde(ms12),'linewidth',2,'color','k');
title('Session 2 vs 1');
xlim([.75 3.25]);
xticks(1:3);
xticklabels({'BL','FU','EOS'});
ylabel('Slope')
set(gca,'linewidth',1,'fontsize',12,'fontweight','bold');

subplot(122);
hold on;
% plot(1:3, log10(SPR_ap), 'linewidth',1,'color',[col col col]);
% plot(1:3, mean(log10(SPR_ap),1),'linewidth',2,'color','k');
plot(1:3, ms13, 'linewidth',1,'color',[col col col]);
errorbar(1:3, mean(ms13,1),withinstde(ms13),'linewidth',2,'color','k');
title('Session 3 vs 1');
xlim([.75 3.25]);
xticks(1:3);
xticklabels({'BL','FU','EOS'});
ylabel('Slope')
set(gca,'linewidth',1,'fontsize',12,'fontweight','bold');


%%
ROI = [12 19 20 48 32 31 49 56 57];

slope_ROI = mean(slope(:,:,ROI),3);

[~,p12_ROI,~,stat12] = ttest(slope_ROI(:,1),slope_ROI(:,2));
[~,p13_ROI,~,stat13] = ttest(slope_ROI(:,1),slope_ROI(:,3));
[~,p23_ROI,~,stat23] = ttest(slope_ROI(:,2),slope_ROI(:,3));

%%
load('fooof_results_chan.mat');
load('ap_removed_chan.mat');
load('fooof_chans_statistics.mat');
load('chans.mat');
figure
subplot(121);
topoplot(squeeze(mean(slope(:,2,:) - slope(:,1,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off', ...
    'emarker2',{find(chans_neg12), '.', 'w', 15},'whitebk','on');
caxis([-.4 .4]);
h = colorbar;
title('Baseline vs Follow up');
ylabel(h,'Slope difference');
set(h,'fontsize',14,'fontweight','bold','linewidth',1.5)
set(gca,'fontsize',14,'fontweight','bold','linewidth',1.5)

subplot(122);
topoplot(squeeze(mean(slope(:,3,:) - slope(:,1,:),1)),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off', ...
    'emarker2',{find(chans_neg13), '.', 'w', 15},'whitebk','on');
caxis([-.4 .4]);
h = colorbar;
title('Baseline vs End of study');
ylabel(h,'Slope difference');
set(h,'fontsize',14,'fontweight','bold','linewidth',1.5);
set(gca,'fontsize',14,'fontweight','bold','linewidth',1.5);

col = .75;
ms13 = mean(slope(:,:,find(chans_neg13)),3);
ms12 = mean(slope(:,:,find(chans_neg12)),3);

% subplot(133);
figure;
hold on;
% plot(1:3, ms13, 'linewidth',1,'color',[col col col]);
errorbar(1:3, mean(ms13,1),withinstde(ms13),'linewidth',2,'color',[.2 .2 .2]);
xlim([.75 3.25]);
xticks(1:3);
xticklabels({'Baseline','Follow up','End of study'});
ylabel('Aperiodic Slope')
set(gca,'linewidth',1.5,'fontsize',15,'fontweight','bold');