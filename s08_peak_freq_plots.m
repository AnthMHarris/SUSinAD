
num_sessions = 3;
sess_ids = {'BL','FU','EOS'};
num_ps = 12;

frex = 3:.1:40;

load('fooof_results.mat');
alpha_lims = [6 12];

peak_alpha = nan(num_ps,num_sessions);

for p = 1:num_ps
    for s = 1:num_sessions
        dat = fooof_results(p,s);
        
        if ~isempty(dat.peak_params) && any(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2))
            if sum(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2))==1
                peak_alpha(p,s) = dat.peak_params(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2),1);
            else
                disp(['Warning: Multifrex p = ' num2str(p) ' s = ' num2str(s)])
                peak_alpha(p,s) = mean(dat.peak_params(dat.peak_params(:,1)>=alpha_lims(1) & dat.peak_params(:,1)<=alpha_lims(2),1));
            end
        end
    end
end

save('peak_alpha.mat',"peak_alpha");

%%
load('peak_alpha.mat');

f = figure; f.Position = [10 10 500 500];
col = .75;
hold on;
plot(1:3, peak_alpha, 'linewidth',1,'color',[col col col]);
errorbar(1:3, mean(peak_alpha, 1, 'omitnan'),withinstde(peak_alpha),'linewidth',2,'color','k');
title('Alpha Peak Frequency - eyes open');
xlim([.75 3.25]);
ylim(alpha_lims);
xticks(1:3);
xticklabels({'BL','FU','EOS'});
ylabel('Peak Alpha Frequency (Hz)')
set(gca,'linewidth',1,'fontsize',12,'fontweight','bold');

%%
[~,p12] = ttest(peak_alpha(:,1), peak_alpha(:,2));
[~,p23] = ttest(peak_alpha(:,2), peak_alpha(:,3));
[~,p13] = ttest(peak_alpha(:,1), peak_alpha(:,3));