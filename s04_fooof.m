
num_sessions = 3;
sess_ids = {'BL','FU','EOS'};
num_ps = 12;

min_frex = 3;
max_frex = 40;

frex = min_frex:.1:max_frex;

load('mean_spect.mat');

py.importlib.import_module('os');
py.importlib.import_module('numpy');
py.importlib.import_module('fooof');

settings = struct(); 
settings.peak_width_limits = [.5 9];
settings.max_n_peaks = 6;
settings.min_peak_height = 0;
settings.peak_threshold = 2.0;
settings.aperiodic_mode = 'fixed'; %'fixed'
settings.verbose = true;
f_range = [min_frex max_frex];

slope = nan(num_ps,num_sessions);
offset = slope;

for p = 1:num_ps

    for s = 1:num_sessions
        session_id = sess_ids{s};

        fooof_results(p,s) = fooof(frex, squeeze(mean_spect(p,s,:)), f_range, settings, true); %#ok<SAGROW> 
        ap_removed(p,s,:) = fooof_results(p,s).power_spectrum - fooof_results(p,s).ap_fit; %#ok<SAGROW> 

        slope(p,s) = fooof_results(p,s).aperiodic_params(2);
        offset(p,s) = fooof_results(p,s).aperiodic_params(1);

    end
end

% %save fooof params and residual spectra separately
save('fooof_results.mat','fooof_results', 'slope', 'offset');
save('ap_removed.mat','ap_removed','frex');

plotfrex = fooof_results(1,1).freqs;
%%
figure; plot(plotfrex, squeeze(mean(ap_removed,1)),'linewidth',1.5);
legend('BL','FU','EOS');

figure; subplot(121); UnivarScatter(slope); subplot(122); UnivarScatter(offset);

[~,p12] = ttest(slope(:,1),slope(:,2));
[~,p13] = ttest(slope(:,1),slope(:,3));
[~,p23] = ttest(slope(:,2),slope(:,3));

