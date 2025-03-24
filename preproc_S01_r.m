%% paths etc.
ftpath     = '/loctmp/AG_Greenlee/m-lib/fieldtrip-20160128';
rootpath   = '/media/Gregor/Experimente/CFS_EEG/';
rawpath    = [rootpath, 'raw/'];
cleanpath  = [rootpath, 'clean/'];
orgpath    = [rootpath, 'org/'];

% filenames
rawname   = 'CFS_EEG0001.vhdr';
rtname    = '1cfsEEG_02-Feb-2017.mat';
vpcode    = 'S01';
dataset   = [rawpath, rawname];
datasetRT = [rawpath, rtname];

% set paths
addpath(ftpath, orgpath); ft_defaults;

%% define segments
cfg=[];
cfg.trialfun = 'trialfun_cfs_responselocked';
cfg.trialdef.prestim    = 2.5; 
cfg.trialdef.poststim   = 1.5;  
cfg.trialdef.eventvalue = 'S  1';
cfg.dataset   = dataset;
cfg.datasetRT = datasetRT;
cfg= ft_definetrial(cfg);

%% preproc eyes
cfgeye = cfg;
cfgeye.demean = 'yes';
cfgeye.detrend = 'yes';
cfgeye.lpfilter = 'yes';
cfgeye.lpfreq = 15;
cfgeye.channel = {'VEOG', 'HEOG'}
eyechans = ft_preprocessing(cfgeye);

cfgdb= [];
cfgdb.viewmode='vertical';
cfgdb.continuous = 'no';
arti1=ft_databrowser(cfgdb, eyechans); 
horizEyeMovements = arti1.artfctdef; clear arti1
fprintf('%i horizontal eye movements are marked.\n', size(horizEyeMovements.visual.artifact,1));
clear eyechans cfgeye cfgdb

%% preproc I
cfg.demean = 'yes';
cfg.detrend = 'yes';
cfg.channel = {'all', '-VEOG', '-HEOG'}
preproc = ft_preprocessing(cfg);

%% channelrepair
badchannels = {'C53', 'C62', 'C3'}; % if necessary; otherwise leave empty 
load ('63equidistant_elec_GV.mat'); % in org-path
cfg=[]
cfg.method = 'distance';
cfg.neighbourdist = 4.9;
cfg.elec = elec;
neighbours = ft_prepare_neighbours(cfg);

cfg=[];
cfg.elec = elec;
cfg.neighbours = neighbours;
cfg.badchannel = badchannels;
interpolated = ft_channelrepair(cfg, preproc); clear preproc


%% clean
cfg= [];
cfg.viewmode='vertical';
cfg.continuous = 'no';
cfg.artfctdef = horizEyeMovements;
cfg = ft_databrowser(cfg, interpolated); 

clean = ft_rejectartifact(cfg, interpolated); clear interpolated

% save
save([cleanpath, 'clean_r_', vpcode], 'clean');
rmpath(orgpath);

