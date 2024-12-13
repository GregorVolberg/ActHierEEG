% randomization for cluster size correction
% H0: true difference between rhos is 0
% test cluster size for random permutations of RDM vector
addpath('./func');
dapath      = './DA/';

%% one-line functions
r2z = @(x) atanh(x); % r to fisher's z
z2r = @(x) tanh(x);  % fisher's z to r

%% load file 
RDMs = get_modelRDMs();
da   = get_decodingaccuracy_fromfiles(dapath);
da_grandmean = squeeze(mean(da.data, 1)); % mean decoding accuracy

nconditions  = size(da.data,2);
ind          = find(tril(ones(nconditions, nconditions), -1) == 1); % lower triangle indices

alpha=0.05;
n = length(ind);
spearmanLowerCI2D = @(x, alpha) reshape(z2r(r2z(x(:)) + norminv(1-(alpha/2)) .* -1 ./ sqrt(1.06*(n-3))), size(x)); % https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation incl. comments
spearmanHigherCI2D = @(x, alpha) reshape(z2r(r2z(x(:)) + norminv(1-(alpha/2)) ./ sqrt(1.06*(n-3))), size(x)); % https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation incl. comments

rng(16);
numperm = 1000;
mxc = {};
for rdm = 1:numel(RDMs)
rdmvector = RDMs{1};
rdmvector = rdmvector(ind);
binarymap = NaN(2, 27,90); 
for permrun = 1:numperm
    rperm = rdmvector(randperm(length(ind)));
    rho = NaN(27,90); 
    for freqbin  = 1:size(da_grandmean,3)
        for tbin = 1:size(da_grandmean, 4)
        da_rdm   = squeeze(da_grandmean(:,:,freqbin, tbin));
        rho(freqbin, tbin) = corr(rperm, da_rdm(ind), 'type', 'Spearman');
        end
    end
binarymap(1,:,:) = rho > 0 & spearmanLowerCI2D(rho, alpha) > 0; %positive
binarymap(2,:,:) = rho < 0 & spearmanHigherCI2D(rho, alpha) < 0; % negatiive
for posneg = 1:2
k = bwlabel(squeeze(binarymap(posneg,:,:)), 4); %4-connected; no vertical
l = regionprops(k);
if isempty(l)
mxc{posneg, rdm, permrun} = 0;
else
mxc{posneg, rdm, permrun} = [l.Area];    
end
end
end
end
save([dapath, 'mxc3rdms.mat'], 'mxc');
exit