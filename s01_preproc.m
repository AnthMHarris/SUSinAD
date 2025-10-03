folders = dir('V:\EEG\*');
folders = folders(3:end);

num_sessions = 3;
sess_ids = {'BL','FU','EOS'};
num_ps = length(folders);

for p = 1:num_ps
    for s = 1:num_sessions
        session_id = sess_ids{s};
        location = [folders(p).folder '\' folders(p).name '\' folders(p).name '_' session_id '\' folders(p).name '_' session_id];
        EEG = pop_readbdf([location '.bdf'], [] ,73,[71 72] );
        EEG=pop_chanedit(EEG, 'lookup','C:\\Users\\anthony\\Documents\\MATLAB\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
        EEG = pop_resample( EEG, 256);
        EEG = pop_eegfiltnew(EEG, 'locutoff',0.3,'plotfreqz',0);
        EEG = pop_eegfiltnew(EEG, 'hicutoff',100,'plotfreqz',0);
        EEG = pop_epoch( EEG, {  '100'  }, [0  300], 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
        dat = cat(3,reshape(EEG.data(:,:,1),72,768,100),reshape(EEG.data(:,:,2),72,768,100),reshape(EEG.data(:,:,3),72,768,100));
        EEG.data = dat;
        EEG.pnts = 768;
        EEG.trials = 300;
        EEG.times = EEG.times(1:768);
        % create sub-epochs
        for t = 1:300
            EEG.epoch(t).event = t;
            EEG.epoch(t).eventlatency = 0;
            EEG.epoch(t).eventtype = 100;
            EEG.epoch(t).eventurevent = t*2-1;
        end
        oldEEG = EEG;
        if p==12 && s==1
             EEG = pop_select( EEG, 'rmchannel',{'PO4','O2'});
        end
        EEG = pop_rejchan(EEG,'elec',1:64,'threshold',5,'measure','kurt','norm','on')       ;
        EEG = pop_rmbase( EEG, [],[]);
        EEG = pop_reref( EEG, [],'exclude',[63:70] );
        EEG = pop_autorej(EEG,'threshold',1000,'electrodes',1:64,'nogui','on');
        r = rank(EEG.data(:,:,1));
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1,'pca',r);

        save([location '_preproc.mat'],'EEG','oldEEG');
    end
end
