addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');

icapath   = './ica';
icacleanedpath = './cleaned';
lay       = '../Data/63equidistant_GreenleeLab_lay.mat';

tmp = dir([icapath, filesep, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

% apply manually to all subjects
clear data_ica_cleaned
vpnr = 22;
icafile = load(files.name{vpnr});

% plot the components for visual inspection
cfg = [];
cfg.component = 1:20;       % specify the component(s) that should be plotted
cfg.layout    = lay; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, icafile.components)
 
% remove the bad components and backproject the data
cfg = [];
cfg.component = [1];% delete e.g. cfg.component = [2]
data_ica_cleaned = ft_rejectcomponent(cfg, icafile.components);
data_ica_cleaned.badcomponents = cfg.component;

save([icacleanedpath, filesep, files.vps{vpnr}, 'cleaned.mat'], 'data_ica_cleaned');


