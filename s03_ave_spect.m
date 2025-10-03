folders = dir('V:\EEG\*');
folders = folders(3:end);
num_sessions = 3;
sess_ids = {'BL','FU','EOS'};
num_ps = length(folders);

frex = 1:.1:40;
mean_spect = nan(num_ps,num_sessions,length(frex));

for p = 1:num_ps
    for s = 1:num_sessions
        session_id = sess_ids{s};
        location = [folders(p).folder '\' folders(p).name '\' folders(p).name '_' session_id '\' folders(p).name '_' session_id];
        
        load([location '_clean.mat']);
        
        spect = nan(64,length(frex));
        for c = 1:64
            spect(c,:) = trimmean(pmtm(squeeze(double(EEG.data(c,:,:))),[],frex,EEG.srate),20,'round',2);

        end

        mean_spect(p,s,:) = mean(spect,1);

        save([location '_spect.mat'],'spect');
    end
end
save('mean_spect.mat','mean_spect','frex');

%%
load('mean_spect.mat');

figure; plot(frex, log10(squeeze(mean(mean_spect,1))),'linewidth',1.5);
legend('BL','FU','EOS');
xlabel('Frequency (Hz)');
ylabel('Log10 Power');
set(gca, 'linewidth',1,'fontsize',12,'fontweight','bold','box','off');

%%
figure; loglog(frex, squeeze(mean(mean_spect,1)),'linewidth',1.5);
legend('BL','FU','EOS'); legend boxoff;
xlabel('log10 Frequency (Hz)');
ylabel('Log10 Power');
set(gca, 'linewidth',1,'fontsize',12,'fontweight','bold','box','off');
xlim([.9 45]);
