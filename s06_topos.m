folders = dir('V:\EEG\*');
folders = folders(3:end);
num_ps = length(folders);
num_sessions = 3;
sess_ids = {'BL','FU','EOS'};

frex = 3:.5:40;

delta = [3 4]; delta_idx = frex >= delta(1) & frex < delta(2);
theta = [4 6]; theta_idx = frex >= theta(1) & frex < theta(2);
alpha = [6 12]; alpha_idx = frex >= alpha(1) & frex < alpha(2);
beta = [12 30]; beta_idx = frex >= beta(1) & frex < beta(2);

load('mean_spect.mat');

delta_power = nan(num_ps,num_sessions,64);
theta_power = delta_power;
alpha_power = delta_power;
beta_power = delta_power;

for p = 1:num_ps
    for s = 1:num_sessions
        session_id = sess_ids{s};
        location = [folders(p).folder '\' folders(p).name '\' folders(p).name '_' session_id '\' folders(p).name '_' session_id];
        load([location '_spect.mat'],'spect');

        delta_power(p,s,:) = mean(spect(:,delta_idx),2);
        theta_power(p,s,:) = mean(spect(:,theta_idx),2);
        alpha_power(p,s,:) = mean(spect(:,alpha_idx),2);
        beta_power(p,s,:) = mean(spect(:,beta_idx),2);
        
    end
end

load('chans.mat');
chans = chans(1:64);

%%


figure; 
for s = 1:num_sessions
    if s==1
        subplot(3,4,((s-1)*4)+1);
        topoplot(squeeze(mean(delta_power(:,s,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        caxis([1 5])
        title('delta');
        subplot(3,4,((s-1)*4)+2);
        topoplot(squeeze(mean(theta_power(:,s,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        caxis([1 5])
        title('theta');
        subplot(3,4,((s-1)*4)+3);
        topoplot(squeeze(mean(alpha_power(:,s,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        title('alpha');
        caxis([0.5 2.5])
        subplot(3,4,((s-1)*4)+4);
        topoplot(squeeze(mean(beta_power(:,s,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        title('beta');
        caxis([0 1.5])
    else
        range = 1.5;
        subplot(3,4,((s-1)*4)+1);
        topoplot(squeeze(mean(delta_power(:,s,:),1))-squeeze(mean(delta_power(:,1,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        caxis([-range range]);
        subplot(3,4,((s-1)*4)+2);
        topoplot(squeeze(mean(theta_power(:,s,:),1))-squeeze(mean(theta_power(:,1,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        caxis([-range range]);
        subplot(3,4,((s-1)*4)+3);
        topoplot(squeeze(mean(alpha_power(:,s,:),1))-squeeze(mean(alpha_power(:,1,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        caxis([-range range]);
        subplot(3,4,((s-1)*4)+4);
        topoplot(squeeze(mean(beta_power(:,s,:),1))-squeeze(mean(beta_power(:,1,:),1)),chans,'numcontour',0,'maplimits','maxmin','gridscale',200,'colormap',parula,'conv','on');
        colorbar
        caxis([-range range]);
    end
end