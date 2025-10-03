clear all; close all; clc

folders = dir('V:\EEG\*');
folders = folders(3:end);

num_sessions = 3;
sess_ids = {'BL','FU','EOS'};

% this part is manual
p = 12;
s = 1;

session_id = sess_ids{s};
location = [folders(p).folder '\' folders(p).name '\' folders(p).name '_' session_id '\' folders(p).name '_' session_id];

load([location '_preproc.mat']);

EEG = SASICA(EEG);

%%

EEG = pop_subcomp(EEG,[],0);
EEG = pop_interp(EEG, oldEEG.chanlocs, 'spherical'); %clear oldEEG

save([location '_clean.mat'],'EEG');

disp('done');