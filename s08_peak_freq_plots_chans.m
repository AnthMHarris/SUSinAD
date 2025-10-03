
num_sessions = 3;
sess_ids = {'BL','FU','EOS'};
num_ps = 12;

frex = 3:.1:40;

load('fooof_results_chan.mat');
alpha_lims = [5 15];

peak_alpha = nan(num_ps,num_sessions,64);

for c = 1:64
    for p = 1:num_ps
        for s = 1:num_sessions
            dat = fooof_results(p,s,c);
            
            if ~isempty(dat.peak_params) && any(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2))
                if sum(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2))==1
                    peak_alpha(p,s,c) = dat.peak_params(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2),1);
                else
                    disp(['Warning: Multifrex p = ' num2str(p) ' s = ' num2str(s)])
                    peak_fs = dat.peak_params(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2),1);
                    peak_pows = dat.peak_params(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2),2);
                    peak_alpha(p,s,c) = peak_fs(peak_pows==max(peak_pows));
                end
            end
        end
    end
end

save('peak_alpha_chan.mat',"peak_alpha");

%%
load('peak_alpha_chan.mat');

load('chans.mat');
figure; 
subplot(231);
topoplot(squeeze(mean(peak_alpha(:,1,:),1,'omitnan')),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([7 10]);
h = colorbar;
title('Baseline');
ylabel(h,'Alpha peak frequency (Hz)');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);
subplot(232);
topoplot(squeeze(mean(peak_alpha(:,2,:),1,'omitnan')),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([7 10]);
h = colorbar;
title('Follow up');
ylabel(h,'Alpha peak frequency (Hz)');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);
subplot(233);
topoplot(squeeze(mean(peak_alpha(:,3,:),1,'omitnan')),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([7 10]);
h = colorbar;
title('End of study');
ylabel(h,'Alpha peak frequency (Hz)');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);

subplot(234);
topoplot(squeeze(mean(peak_alpha(:,2,:) - peak_alpha(:,1,:),1,'omitnan')),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([-1 1]);
h = colorbar;
title('Baseline vs Follow up');
ylabel(h,'Alpha peak difference (Hz)');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);
subplot(235);
topoplot(squeeze(mean(peak_alpha(:,3,:) - peak_alpha(:,1,:),1,'omitnan')),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([-1 1]);
h = colorbar;
title('Baseline vs End of study');
ylabel(h,'Alpha peak difference (Hz)');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);
subplot(236);
topoplot(squeeze(mean(peak_alpha(:,3,:) - peak_alpha(:,2,:),1,'omitnan')),chans(1:64),'numcontour',0,...
    'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','off','whitebk','on');
caxis([-1 1]);
h = colorbar;
title('Follow up vs End of study');
ylabel(h,'Alpha peak difference (Hz)');
set(h,'fontsize',12,'fontweight','bold','linewidth',1);
set(gca,'fontsize',12,'fontweight','bold','linewidth',1);

%%
[~,p12,~,stat12] = ttest(squeeze(peak_alpha(:,1,:)),squeeze(peak_alpha(:,2,:)));
[~,p13,~,stat13] = ttest(squeeze(peak_alpha(:,1,:)),squeeze(peak_alpha(:,3,:)));
[~,p23,~,stat23] = ttest(squeeze(peak_alpha(:,2,:)),squeeze(peak_alpha(:,3,:)));

disp(min(bonf_holm(p12)));
disp(min(bonf_holm(p13)));
disp(min(bonf_holm(p23)));

load('chans.mat');
nei = spatial_neighbors(chans(1:64),56);

[p_pos12, chans_pos12, p_neg12, chans_neg12] = ClustPerm1D(p12', nei, mean(squeeze(peak_alpha(:,1,:)) - squeeze(peak_alpha(:,2,:)),1)');
[p_pos13, chans_pos13, p_neg13, chans_neg13] = ClustPerm1D(p13', nei, mean(squeeze(peak_alpha(:,1,:)) - squeeze(peak_alpha(:,3,:)),1)');
[p_pos23, chans_pos23, p_neg23, chans_neg23] = ClustPerm1D(p23', nei, mean(squeeze(peak_alpha(:,2,:)) - squeeze(peak_alpha(:,3,:)),1)');

%%
out_12 = cluster_perm_1d(squeeze(peak_alpha(:,2,:)-peak_alpha(:,1,:)),nei);
out_13 = cluster_perm_1d(squeeze(peak_alpha(:,3,:)-peak_alpha(:,1,:)),nei);
out_23 = cluster_perm_1d(squeeze(peak_alpha(:,3,:)-peak_alpha(:,2,:)),nei);


%%
focal_chans = [12 19 20 48 32 31 49 56 57];


focal_alpha = mean(peak_alpha(:,:,focal_chans),3,'omitnan');

[~,focal_p12] = ttest(focal_alpha(:,1), focal_alpha(:,2));
[~,focal_p13] = ttest(focal_alpha(:,1), focal_alpha(:,3));
[~,focal_p23] = ttest(focal_alpha(:,2), focal_alpha(:,3));

f = figure; f.Position = [10 10 500 500];
col = .75;
hold on;
plot(1:3, focal_alpha, 'linewidth',1,'color',[col col col]);
errorbar(1:3, mean(focal_alpha, 1, 'omitnan'),withinstde(focal_alpha),'linewidth',2,'color','k');
title('Alpha Peak Frequency - eyes open');
xlim([.75 3.25]);
ylim(alpha_lims);
xticks(1:3);
xticklabels({'BL','FU','EOS'});
ylabel('Peak Alpha Frequency (Hz)')
set(gca,'linewidth',1,'fontsize',12,'fontweight','bold');