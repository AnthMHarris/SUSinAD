folders = dir('V:\EEG\*');
folders = folders(3:end);
num_ps = length(folders);
num_sessions = 3;
sess_ids = {'BL','FU','EOS'};

frex = [8 14];
mean_cohere = nan(num_ps,num_sessions);

times2save = [500 2500];
        

for p = 1:num_ps
    for s = 1:num_sessions
        session_id = sess_ids{s};
        location = [folders(p).folder '\' folders(p).name '\' folders(p).name '_' session_id '\' folders(p).name '_' session_id];
        load([location '_clean.mat']);

        times2saveidx = EEG.times >= times2save(1) & EEG.times <= times2save(2);

        %filter parameters
        nyquist = EEG.srate/2;
        lower_filter_bound = 6; % Hz 8
        upper_filter_bound = 12; % Hz 14
        transition_width   = 0.25;
        filter_order       = round(3*(EEG.srate/lower_filter_bound));
        
        %make filter
        ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
        idealresponse = [ 0 0 1 1 0 0 ];
        filterweights = firls(filter_order,ffrequencies,idealresponse);
        
        complex_fd = nan(64,768,size(EEG.data,3));
        for t = 1:size(EEG.data,3)
            for c = 1:64
                filtered_data = filtfilt(filterweights,1,double(EEG.data(c,:,t)));
                complex_fd(c,:,t) = hilbert(filtered_data);
            end
        end
        %compute pairwise coherence
        coh = nan(64,64);
        for c1 = 1:63
            for c2 = c1+1:64
                sig1 = squeeze(complex_fd(c1,:,:));
                sig2 = squeeze(complex_fd(c2,:,:));

                % compute power and cross-spectral power
                spec1 = mean(sig1.*conj(sig1),2);
                spec2 = mean(sig2.*conj(sig2),2);
                specX = abs(mean(sig1.*conj(sig2),2)).^2;


                % compute spectral coherence, using only requested time points
                coh(c1,c2) = mean(specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx)));

            end
        end
    

        mean_cohere(p,s) = mean(coh,'all','omitnan');

    end
end

save('mean_cohere.mat',"mean_cohere");

%%
load('mean_cohere.mat');

[~,p12] = ttest(mean_cohere(:,1),mean_cohere(:,2));
[~,p13] = ttest(mean_cohere(:,1),mean_cohere(:,3));
[~,p23] = ttest(mean_cohere(:,2),mean_cohere(:,3));


f = figure; f.Position = [10 10 500 500];
col = .75; 
hold on;
plot(1:3, mean_cohere, 'linewidth',1,'color',[col col col]);
errorbar(1:3, mean(mean_cohere),withinstde(mean_cohere),'linewidth',2,'color','k');
title('Alpha Coherence');
xlim([.75 3.25]);
xticks(1:3);
xticklabels({'BL','FU','EOS'});
ylabel('Spectral Coherence')
set(gca,'linewidth',1,'fontsize',12,'fontweight','bold');