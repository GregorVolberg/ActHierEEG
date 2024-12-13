% addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');
%addpath(genpath('./func/libsvm'));

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

erppath    = './erp/';
tmp = dir([erppath, '*erp.mat']);
tmpchar = char({tmp.name});
tfrfiles.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
tfrfiles.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

dapath     = './DAerp/';
tmp = dir([dapath, '*DAmean.mat']);
tmpchar = char({tmp.name});
dafiles.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
dafiles.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

%% subject loop
% for TFR
% baseline  = [-0.5 -0.1];
% dummy     = load(tfrfiles.name{1});
% bsl       = dsearchn(dummy.tfr.time', baseline');
% bslInd    = bsl(1):bsl(2);
% taxisTFR  = dummy.tfr.time;
% clear dummy
% for vp = 1:numel(tfrfiles.name)
% tmp       = load (tfrfiles.name{vp});
% tmpm      = squeeze(mean(tmp.tfr.powspctrm, 1));
% bslpow    = nanmean(tmpm(:, :, bslInd),3); % 62 x 27
% dbpow(vp,:,:,:) = 10*log10(tmpm ./ bslpow);
% clear tmp
% end
% meandbpow = squeeze(mean(mean(dbpow, 1), 2)); % across channels and conditions
% save([erppath, 'meandbpow.mat'], 'meandbpow');

% 3D-DA with single subjects
for vp = 1:numel(dafiles.name)
da = load (dafiles.name{vp});
DA_all(vp, :, :) = squeeze(nanmean(nanmean(da.DAmean,1), 2));
end
save([dapath, 'DA_all.mat'], 'DA_all');

%% stats: randomization test per bin, one sample versus 50
% step 1: determine p-value per bin
% - create permutation t distribution
% - apply to t-values of actual data
rng(22);
numruns = 10000;
tmpt = nan(numruns,600);
tmpp = nan(numruns,600);
for permruns = 1:numruns
permVector = get_permvector(size(DA_all,1)); 
tmpdam = abs(DA_all-permVector);
[~, p, ~, t] = ttest(tmpdam, 50);
tmpt(permruns, :) = squeeze(t.tstat);
tmpp(permruns, :) = squeeze(p);
end
[~, p, ~, t] = ttest(DA_all, 50);
pnum = squeeze(sum((tmpt > 0) & (tmpt > t.tstat),1)); 
pnum(isnan(squeeze(p))) = nan; % exclude nan bins at 4 Hz
pval = pnum / numruns;


permttest.numruns    = numruns;   
permttest.tmat       = tmpt;
permttest.pmat       = tmpp;
permttest.t          = t; % actual t value
permttest.p          = p; % analytical p value
permttest.pval       = pval; % permutation p value
save([dapath, 'permttest.mat'], 'permttest');

%% stats step 2: thresholding
% - bonferroni correction 
% - check for number of connected bins (with positive t)
% use only post-stimulus interval (0:1.18s)
pthresh = 0.025 %/ (27*90); % bonferroni .05, on all bins
bonfThreshedP = (pval) < pthresh;
[actualP, N] = bwlabel(bonfThreshedP, 4); % 4-connected (no vertical connection)
clusters    = arrayfun(@(x) find(actualP==x), 1:N, 'UniformOutput', false); % sorted labels 1:4
clusterSize = arrayfun(@(x) numel([x{:}]), clusters);

%HzRange = [1:27]; % 4:30 Hz
mxc = cell(1, size(tmpp,1));
tmppPos = tmpp;
tmppPos(tmpt < 0) = 1; % only use right side t-tests

for nmr = 1:size(tmpp,1)
threshedP = squeeze(tmppPos(nmr, 100:600)) < pthresh; %bonferoni thres
k = bwlabel(threshedP, 4); %4-connected; no vertical
l = regionprops(k);
if isempty(l)
mxc{nmr} = 0;
else
mxc{nmr} = [l.Area];
end
end

% find clusters with p < .05
pSigClusters = arrayfun(@(x) sum([mxc{:}] > x) / numel([mxc{:}]), clusterSize);
sigIndex = find(pSigClusters < .05);

%  properties of significant clusters
tvec = -0.2:0.002:1;
%fvec = 4:1:30;
tmap = squeeze(t.tstat);
damap = squeeze(mean(DA_all,1));
clustermat = [];
for clst = 1:numel(sigIndex)
[r, c] = ind2sub(size(actualP), clusters{sigIndex(clst)});
indx   = [clusters{sigIndex(clst)}];
[maxt, maxtind]   = max(tmap(indx));
[maxda, maxdaind] = max(damap(indx));
clustermat(clst).cluster     = [r, c, indx, tmap(indx), damap(indx)];
clustermat(clst).clustersize = length(r);
clustermat(clst).cluster_p   = pSigClusters(sigIndex(clst));
clustermat(clst).time        = tvec([min(c), max(c)]);
%clustermat(clst).freq        = fvec([min(r), max(r)]);
clustermat(clst).max_t       = [tvec(c(maxtind)), maxt];
clustermat(clst).max_da      = [tvec(c(maxdaind)), maxda];
fprintf('\nCluster %i, %i bins, p = %.4f', clst, length(r), pSigClusters(sigIndex(clst)));
fprintf('\n %.4f to %.4fs', tvec(min(c)), tvec(max(c)));
%fprintf('\n %i to %.i Hz', fvec(min(r)), fvec(max(r)));
fprintf('\n max t  = %.3f at %.4fs', maxt, tvec(c(maxtind)));
fprintf('\n max DA = %.3f at %.4fs\n\n', maxda, tvec(c(maxdaind)));
end

save([dapath, 'clustermat.mat'], 'clustermat');

%% Results

% Cluster 1, 30 bins, p = 0.0460
%  0.0000 to 0.4600s
%  4 to 5 Hz
%  max t  = 10.646 at 0.3000s and 5 Hz
%  max DA = 56.543 at 0.2800s and 5 Hz
% 
% 
% Cluster 2, 27 bins, p = 0.0496
%  0.0800 to 0.3400s
%  7 to 9 Hz
%  max t  = 9.869 at 0.2200s and 7 Hz
%  max DA = 58.129 at 0.2000s and 7 Hz
