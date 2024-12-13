addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');

logpath   = '../Data/log/';
rawpath   = '../Data/raw/';
cleanpath = '../Results/Preprocessing/';
lay       = '../Data/63equidistant_GreenleeLab_lay.mat';

% get files (in ./func)
ffiles = getfiles(logpath, rawpath, cleanpath);

% prepare layout 
cfg=[]; 
cfg.layout = lay;
layout = ft_prepare_layout(cfg);

for vpnr = 1:numel(ffiles.vhdr)
%% trialdef
cfg                         = [];
cfg.dataset                 = ffiles.vhdr{vpnr};
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'Stimulus';
cfg.trialdef.eventvalue     = {'S  1', 'S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12'};
cfg.trialdef.prestim        = 1.2; % in seconds 
cfg.trialdef.poststim       = 1.8; % in seconds
cfg = ft_definetrial(cfg);

% subjects [1:4, 6:10] have one false extra marker in condition 1 
% directly at the beginning of the experiment (and so 865 instead of 864
% trials). Remove first row of trial definition for these cases
targetSubjects = arrayfun(@(x) strcat('SUB', num2str(x,'%02.f')), [1:4, 6:10], 'UniformOutput', false);
if ismember(ffiles.vps{vpnr}, targetSubjects)
   cfg.trl(1,:) = [];
end

% add trialinfo from protocol file and check consistency with EEG marker values
infomat = load(ffiles.info{vpnr});
tcode = arrayfun(@(x) infomat.ExpInfo.TrialInfo(x).trial.code, 1:numel(infomat.ExpInfo.TrialInfo));
tisi  = arrayfun(@(x) infomat.ExpInfo.TrialInfo(x).trial.pageDuration(2) * 1/60, 1:numel(infomat.ExpInfo.TrialInfo));
tpic  = arrayfun(@(x) infomat.ExpInfo.TrialInfo(x).trial.pageNumber(3), 1:numel(infomat.ExpInfo.TrialInfo));
tindices = find(tcode >= 1 & tcode <= 12);
cfg.trl = [cfg.trl, tcode(tindices)', tisi(tindices)', tpic(tindices)'];
 if any(cfg.trl(:,4) ~= cfg.trl(:,5))
     error('Inconsistent marker values and protocol file\n'); end

% find trials in cfg that correspond to trials in cleaned timelock data and
% select clean trials
cleanTL = load(ffiles.clean{vpnr});
ttrials = arrayfun(@(x) find(x > cfg.trl(:, 1) & x < cfg.trl(:, 2)), cleanTL.timelock.sampleinfo(:,1));
cfg.trl = cfg.trl(ttrials,:);

% add proprocessing options and preprocess
cfg.channel = {'all', '-VEOG', '-HEOG'};
cfg.demean = 'yes'; 
cfg.reref = 'yes'; 
cfg.refchannel    = {'all'};
cfg.hpfilter = 'yes';
cfg.hpfreq = 0.1;
cfg.hpfiltord=5;
pre_proc = ft_preprocessing(cfg);

% % ica
cfg = [];
cfg.channel = 'all';
cfg.method = 'runica';
cfg.numcomponent = 'all';
components = ft_componentanalysis(cfg, pre_proc);
save(['.', filesep, 'ica', filesep, ffiles.vps{vpnr}, 'ica.mat'], 'components');


end

% 864 or 865 trials;
% 12 conditions, 72 trials per condition
% 144 trials with code 13 (catch trials)
% 5 trials with code 14 (?)


