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

rng(16);
numperm = 1000;
[rdm1vector, rdm2vector] = deal(RDMs{1}, RDMs{2});
[rdm1vector, rdm2vector] = deal(rdm1vector(ind), rdm2vector(ind));
rho1 = NaN; rho2 = NaN; mxc = {};
for permrun = 1:numperm
    rndorder = randperm(length(ind));
    [r1, r2] = deal(rdm1vector(rndorder), rdm2vector(rndorder));
    diffmap = nan(size(da_grandmean,3), size(da_grandmean,4));
    for freqbin  = 1:size(da_grandmean,3)
        for tbin = 1:size(da_grandmean, 4)
        da_rdm   = squeeze(da_grandmean(:,:,freqbin, tbin));
        rho1  = corr(r1, da_rdm(ind), 'type', 'Spearman');
        rho2  = corr(r2, da_rdm(ind), 'type', 'Spearman');
 
        diffmap(freqbin, tbin) = z2r(r2z(rho1)*-1 + r2z(rho2)*1);
        end
    end
pzmap = rho2z(diffmap);
binarymap = pzmap >= 1.96;
k = bwlabel(binarymap, 4); %4-connected; no vertical
l = regionprops(k);
if isempty(l)
mxc{permrun} = 0;
else
mxc{permrun} = [l.Area];
end
end
save([dapath, 'mxczmap.mat'], 'mxc');
