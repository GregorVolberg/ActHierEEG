addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');

% paths and file names
cleanpath  = './cleaned/';
itcpath    = './itc/';

tmp = dir([cleanpath, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

% tfr
cfgtfr = [];
cfgtfr.output             = 'fourier';
cfgtfr.method             = 'mtmconvol';
cfgtfr.taper              = 'hanning';
cfgtfr.foi                = 4:1:30; % 4 to 40 Hz
cfgtfr.t_ftimwin          = 5./cfgtfr.foi;
cfgtfr.toi                = -0.6:0.01:1.2;%
cfgtfr.pad                = 'nextpow2';
cfgtfr.keeptrials         = 'yes';

%cfgbsl = [];
%cfgbsl.baseline           = [-0.5 -0.1];
%cfgbsl.baselinetype       = 'db';

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
eeg = load (files.name{vp});
itc =  ft_freqanalysis(cfgtfr, eeg.data_ica_cleaned);
save([itcpath, files.vps{vp}, 'itc.mat'], 'itc', '-v7.3');
clear eeg itc
end % end vp

% see https://www.sciencedirect.com/science/article/pii/S0960982222017729?via%3Dihub
% downsample to 0.02 s bins and baseline correct to db
