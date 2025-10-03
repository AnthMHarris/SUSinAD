%%

load('fooof_results_chan.mat');

frex = 3:.1:40;

delta = [3 4]; delta_idx = frex >= delta(1) & frex < delta(2);
theta = [4 8]; theta_idx = frex >= theta(1) & frex < theta(2);
alpha = [8 13]; alpha_idx = frex >= alpha(1) & frex < alpha(2);
beta = [13 30]; beta_idx = frex >= beta(1) & frex < beta(2);

full_spect = nan(12,3,64,371);
ap_corr = full_spect;
for p = 1:12
    for s = 1:3
        for c = 1:64
            full_spect(p,s,c,:) = fooof_results(p,s,c).power_spectrum;
            ap_corr(p,s,c,:) = fooof_results(p,s,c).power_spectrum - fooof_results(p,s,c).ap_fit;
        end
    end
end

delta_power = mean(full_spect(:,:,:,delta_idx),4);
theta_power = mean(full_spect(:,:,:,theta_idx),4);
alpha_power = mean(full_spect(:,:,:,alpha_idx),4);
beta_power = mean(full_spect(:,:,:,beta_idx),4);

delta_power_ap = mean(ap_corr(:,:,:,delta_idx),4);
theta_power_ap = mean(ap_corr(:,:,:,theta_idx),4);
alpha_power_ap = mean(ap_corr(:,:,:,alpha_idx),4);
beta_power_ap = mean(ap_corr(:,:,:,beta_idx),4);

SPR = (alpha_power + beta_power) ./ (delta_power + theta_power);
SPR_ap = (alpha_power_ap + beta_power_ap) ./ (delta_power_ap + theta_power_ap);


%%
[~,p12] = ttest(squeeze(SPR(:,1,:)), squeeze(SPR(:,2,:))); 
[~,p13] = ttest(squeeze(SPR(:,1,:)), squeeze(SPR(:,3,:))); 
[~,p23] = ttest(squeeze(SPR(:,2,:)), squeeze(SPR(:,3,:))); 

[~,p12_ap] = ttest(squeeze(SPR_ap(:,1,:)), squeeze(SPR(:,2,:))); 
[~,p13_ap] = ttest(squeeze(SPR_ap(:,1,:)), squeeze(SPR(:,3,:))); 
[~,p23_ap] = ttest(squeeze(SPR_ap(:,2,:)), squeeze(SPR(:,3,:))); 

load('chans.mat');
old_path = addpath('C:\Users\uqaharr7\Desktop\mass_univ03272017');
nei = spatial_neighbors(chans(1:64),56);

[p_pos12, chans_pos12, p_neg12, chans_neg12] = ClustPerm1D(p12', nei, mean(squeeze(SPR(:,1,:)) - squeeze(SPR(:,2,:)),1)');
[p_pos13, chans_pos13, p_neg13, chans_neg13] = ClustPerm1D(p13', nei, mean(squeeze(SPR(:,1,:)) - squeeze(SPR(:,3,:)),1)');
[p_pos23, chans_pos23, p_neg23, chans_neg23] = ClustPerm1D(p23', nei, mean(squeeze(SPR(:,2,:)) - squeeze(SPR(:,3,:)),1)');

[p_pos12_ap, chans_pos12_ap, p_neg12_ap, chans_neg12_ap] = ClustPerm1D(p12_ap', nei, mean(squeeze(SPR_ap(:,1,:)) - squeeze(SPR_ap(:,2,:)),1)');
[p_pos13_ap, chans_pos13_ap, p_neg13_ap, chans_neg13_ap] = ClustPerm1D(p13_ap', nei, mean(squeeze(SPR_ap(:,1,:)) - squeeze(SPR_ap(:,3,:)),1)');
[p_pos23_ap, chans_pos23_ap, p_neg23_ap, chans_neg23_ap] = ClustPerm1D(p23_ap', nei, mean(squeeze(SPR_ap(:,2,:)) - squeeze(SPR_ap(:,3,:)),1)');


%%
ROI = [12 19 20 48 32 31 49 56 57];
SPR_ROI = (mean(alpha_power(:,:,ROI),3) + mean(beta_power(:,:,ROI),3)) ./ (mean(delta_power(:,:,ROI),3) + mean(theta_power(:,:,ROI),3));
SPR_ap_ROI = (mean(alpha_power_ap(:,:,ROI),3) + mean(beta_power_ap(:,:,ROI),3)) ./ (mean(delta_power_ap(:,:,ROI),3) + mean(theta_power_ap(:,:,ROI),3));

[~,p12_ROI] = ttest(SPR_ROI(:,1), SPR_ROI(:,2,:)); 
[~,p13_ROI] = ttest(SPR_ROI(:,1,:), SPR_ROI(:,3,:)); 
[~,p23_ROI] = ttest(SPR_ROI(:,2,:), SPR_ROI(:,3,:)); 

[~,p12_ap_ROI] = ttest(SPR_ap_ROI(:,1,:), SPR_ROI(:,2,:)); 
[~,p13_ap_ROI] = ttest(SPR_ap_ROI(:,1,:), SPR_ROI(:,3,:)); 
[~,p23_ap_ROI] = ttest(SPR_ap_ROI(:,2,:), SPR_ROI(:,3,:)); 
