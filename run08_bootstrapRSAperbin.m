% addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;
% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

addpath('./func');
dapath      = './DA/';

% get model RDMS and decoding accuracy
RDMs = get_modelRDMs();
da   = get_decodingaccuracy_fromfiles(dapath);
da_grandmean = squeeze(mean(da.data, 1));

% indices of lower triangle of RDMs
nconditions = size(da_grandmean, 1);
ind = find(tril(ones(nconditions, nconditions), -1) == 1);

rho   = nan(numel(RDMs), size(da_grandmean, 3), size(da_grandmean, 4));
for num_RDM = 1:numel(RDMs)
    RDM = RDMs{num_RDM};
    for freqbin = 1:size(da_grandmean,3)
        for tbin = 1:size(da_grandmean, 4)
        da_rdm = squeeze(da_grandmean(:,:,freqbin, tbin));
        rho(num_RDM, freqbin, tbin) = corr(RDM(ind), ...
                                           da_rdm(ind), ...
                                           'type', 'Spearman');    
        end
    end
end

% bootstrap CI
n_randomizations = 1000;
rho_boot   = nan(n_randomizations, numel(RDMs), size(da_grandmean, 3), size(da_grandmean, 4));
rng(11); % initialize random number
for boot = 1:1000
rnd    = randsample(numel(da.subjects), numel(da.subjects), 'true');
boot_da = squeeze(mean(da.data(rnd, :,:,:,:),1));
for num_RDM = 1:length(RDMs)
    RDM = RDMs{num_RDM};
    for freqbin = 1:size(boot_da,3)
        for tbin = 1:size(boot_da, 4)
        boot_rdm = squeeze(boot_da(:,:,freqbin, tbin));
        rho_boot(boot, num_RDM, freqbin, tbin) = corr(RDM(ind), boot_rdm(ind), 'type', 'Spearman');    
        end
    end
end
end

rho_per_bin.rho = rho;
rho_per_bin.rho_boot = rho_boot;
save([dapath,'rho_per_bin.mat'], 'rho_per_bin');

