addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');

% paths and file names
cleanpath  = './cleaned/';
erppath    = './erp/';

tmp = dir([cleanpath, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

% a la Gregor
cfgpp=[];
cfgpp.channel = {'all', '-VEOG', '-HEOG'};
cfgpp.demean = 'yes'; 
cfgpp.reref = 'yes'; 
cfgpp.refchannel    = {'all'};
cfgpp.hpfilter = 'yes';
cfgpp.hpfreq = 0.1;
cfgpp.hpfiltord=5;

cfgtl = [];
cfgtl.keeptrials = 'yes';

cfgrd = [];
cfgrd.toilim = [-0.2 1];

cfgbsl = [];
cfgbsl.baseline = [-0.2 0];
cfgbsl.parameter = 'trial';

% trial selection
% condition matrix in *.trialinfo
% col1:= condition code, taken from Tonghe's segmented EEG data
% col2:= condition code, taken from protocol file (as consistency check)
% col3:= pre-stimulus ISI in s, computed from protocol file. Was randomized 1 - 1.5 s
% col4:= number refering to target stimulus (see protocol files, e. g. ../Data/log/SUB01_01.mat, ExpInfo.stimNames) 
% condition codes   1:4 are basic level actions locomotion
%                   5:8 are basic level actions ingestion
%                   9:12 are basic level actions cleaning

% subject loop
for vp = 1:numel(files.name) % subject loop
eeg      = load(files.name{vp}, 'data_ica_cleaned');
pre_proc = ft_preprocessing(cfgpp, eeg.data_ica_cleaned);
tl       = ft_timelockanalysis(cfgtl, pre_proc);
redef    = ft_redefinetrial(cfgrd, tl);
erp      = ft_timelockbaseline(cfgbsl,redef);
save([erppath, files.vps{vp}, 'erp.mat'], 'erp', '-v7.3');
clear eeg erp
end % end vp

