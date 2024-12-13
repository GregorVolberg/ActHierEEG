% ==============
% get and plot RSA results, per bin
% ==============

%% paths and files
addpath('./func');
dapath      = './DA/';
rhobootstrapfile = [dapath, 'rho_per_bin.mat'];  

%% one-line functions
r2z = @(x) atanh(x); % r to fisher's z
z2r = @(x) tanh(x);  % fisher's z to r

%% load file 
RDMs = get_modelRDMs();
da   = get_decodingaccuracy_fromfiles(dapath);

%% analysis II: RDM correlation per TF bin
da_grandmean = squeeze(mean(da.data, 1)); % on mean decoding accuracy
nconditions  = size(da.data,2);
ind          = find(tril(ones(nconditions, nconditions), -1) == 1); % lower triangle indices
rho          = nan(numel(RDMs), size(da_grandmean, 3), size(da_grandmean, 4));
p_rho        = nan(numel(RDMs), size(da_grandmean, 3), size(da_grandmean, 4));
for num_RDM = 1:numel(RDMs)
    RDM = RDMs{num_RDM};
    for freqbin = 1:size(da_grandmean,3)
        for tbin = 1:size(da_grandmean, 4)
        da_rdm = squeeze(da_grandmean(:,:,freqbin, tbin));
        [tmp_rho, tmp_p] = corr(RDM(ind), ...
                                da_rdm(ind), ...
                                'type', 'Spearman');    
        rho(num_RDM, freqbin, tbin)   = tmp_rho;
        p_rho(num_RDM, freqbin, tbin) = tmp_p;
        end
    end
end

% plot
titles = {'Subordinate Level', 'Basic Level', 'Superordinate Level'};
scales = {[0,0.15], [0,0.15], [0,0.15]};
figure(1);
for j = 1:3
subplot(3,1,j);
imagesc([-0.6:0.02:1.18],30:-1:4, squeeze(rho(j,:,:)), scales{j});
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Spearman''s rho';
title(titles{j});
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
end
rectFig = get(gcf,'position');
width=600;
height=900;
set(gcf,'position',[rectFig(1),rectFig(2)- height,width,height], 'color', 'white');

% find maximum rho per rdm
timeaxis = [-0.6:0.02:1.18];
freqaxis = [4:30];
[m, itime] = arrayfun(@(x) max(max(squeeze(rho(x,:,:)))), 1:3); % time points
[m, ifreq] = arrayfun(@(x) max(max(squeeze(rho(x,:,:))')), 1:3); % frequencies
[m; freqaxis(ifreq); timeaxis(itime)]'

% find ci 
% Bonnet & Wright (2000), Psychometrika, 65(1), 23-28
% rs = 0.7;
% alph = 0.01;
% n = length(ind);
% b = 3;
% L1 = 0.5*(log(1+rs) - log(1-rs)) - (((1+rs^2/2)^(1/2) * norminv(1-(alph/2))) / ((n-b)^(1/2)));
% L2 = 0.5*(log(1+rs) - log(1-rs)) + (((1+rs^2/2)^(1/2) * norminv(1-(alph/2))) / ((n-b)^(1/2)));
% cilow  = (exp(2*L1)-1) / (exp(2*L1)+1);
% cihigh = (exp(2*L2)-1) / (exp(2*L2)+1);
% [cilow, cihigh]
% rs=rho(rdm, :, itime(rdm))
% ci=bonnetwright(rs, alph, n);


%% find bootstrap ci 
% n_bootstrap = 1000;
% rng(12); % initialize random number 
% rnd               = arrayfun(@(x) randsample(size(da.data, 1), size(da.data, 1), 'true'), ...
%                       1:n_bootstrap, 'UniformOutput', false);
% 
% ... per frequency, at time point with maximum rho                  
% for j = 1:3
% RDM = RDMs{j};
% bootstrapped_da   = arrayfun(@(x) squeeze(mean(da.data(rnd{x},:, :,:,itime(j)))), ...
%                       1:n_bootstrap, 'UniformOutput', false);
% vectorized_da     = arrayfun(@(x) reshape(bootstrapped_da{x}, [72*72, 27]), ...
%                       1:n_bootstrap, 'UniformOutput', false);
% lower_triangle    = cellfun(@(x) x(ind,:), vectorized_da, 'UniformOutput', false);
% 
% SpearmansRho{j}   = cell2mat(arrayfun(@(x) corr(lower_triangle{x}, RDM(ind), 'type', 'Spearman'), ...
%                     1:n_bootstrap, 'UniformOutput', false))';
% end
% for export to R
% for CIs see
% https://link.springer.com/article/10.1007/BF02294183#preview
%% plot
% find ci from z score 
alph = 0.01;
zcrit = norminv(1-(alph/2));
plusminus = [-1,1];
spearmanCI = @(x) tanh(atanh(x) + zcrit * plusminus / sqrt(1.06*(n-3))); % https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation incl. comments

for rdm = 1:3
rs = rho(rdm, :, itime(rdm));
ci = spearmanCI(rs');%[quantile(SpearmansRho{rdm}, twosidedp/2); ...
                    %quantile(SpearmansRho{rdm}, 1- twosidedp/2)]';
figure(rdm+1);                
poly = polyshape([[-1*ci(:,1)'; 4:30]'; [-1*flipud(ci(:,2))'; 30:-1:4]']);
plot(poly, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]); hold on
plot(rs*-1, freqaxis, 'LineWidth', 1.5, 'Color', 'b');
title([num2str(timeaxis(itime(rdm))), ' ms']);
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
xlabel('Spearman''s rho'); ylabel('Frequency (Hz)');
ylim([4 30]); xlim([-0.21 0.21])
set(gca,'YTick',[5:5:30], ...
    'YTickLabel', [5:5:30]);
set(gca,'XTick',[-0.2:0.1:0.2], ...
    'XTickLabel', [0.2:-0.1:-0.2]);
rectFig = get(gcf,'position');
width=150;
height=300;
set(gcf,'position',[rectFig(1),rectFig(2)-height/2,width,height], 'color', 'white');
text(repmat(-0.19, sum(ci(:,1) > 0),1), freqaxis(ci(:,1) > 0), '*', 'FontSize', 14, 'VerticalAlignment', 'middle')
end
 
% ... per time point, at frequency with maximum rho                  
% for j = 1:3
% RDM = RDMs{j};
% bootstrapped_da   = arrayfun(@(x) squeeze(mean(da.data(rnd{x},:, :,ifreq(j),:))), ...
%                       1:n_bootstrap, 'UniformOutput', false);
% vectorized_da     = arrayfun(@(x) reshape(bootstrapped_da{x}, [72*72, 90]), ...
%                       1:n_bootstrap, 'UniformOutput', false);
% lower_triangle    = cellfun(@(x) x(ind,:), vectorized_da, 'UniformOutput', false);
% 
% SpearmansRho{j}   = cell2mat(arrayfun(@(x) corr(lower_triangle{x}, RDM(ind), 'type', 'Spearman'), ...
%                     1:n_bootstrap, 'UniformOutput', false))';
% end

% plot
for rdm = 1:3
rs = squeeze(rho(rdm, ifreq(rdm), :));
ci = spearmanCI(rs);%[quantile(SpearmansRho{rdm}, twosidedp/2); ...
                    %quantile(SpearmansRho{rdm}, 1- twosidedp/2)]';

%ci = [quantile(SpearmansRho{rdm}, twosidedp/2); ...
%                    quantile(SpearmansRho{rdm}, 1- twosidedp/2)]';
figure(rdm+1);                
poly = polyshape([[timeaxis; ci(:,1)']'; [fliplr(timeaxis); flipud(ci(:,2))']']);
plot(poly, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]); hold on
plot(timeaxis, rs, 'LineWidth', 1.5, 'Color', 'b');
title([num2str(freqaxis(ifreq(rdm))), ' Hz']);
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
xlabel('Time (s)'); ylabel('Spearman''s rho');
ylim([-0.1 0.21]); xlim([-0.6 1.2])
% set(gca,'YTick',[-0.1:0.1:0.2], ...
%     'YTickLabel', [-0.1:0.1:0.2]);
set(gca,'XTick',[-0.6:0.2:1.2], ...
    'XTickLabel', [-0.6:0.2:1.2]);
% 
rectFig = get(gcf,'position');
width=300;
height=100;
set(gcf,'position',[rectFig(1),rectFig(2)-height/2,width,height], 'color', 'white');
text(timeaxis(ci(:,1) > 0), repmat(0.2, sum(ci(:,1) > 0),1), '*', 'FontSize', 14, 'VerticalAlignment', 'middle')
end


% difference subordinate vs basic (i.e., [-0.5 0.5 0])
diffmap = z2r(squeeze(r2z(rho(1,:,:)))*-1 + squeeze(r2z(rho(2,:,:)))*1);
rho2z = @(x) sqrt((numel(ind)-3)/1.06) * x; % transform Fisher zr into z for obtaining p, see https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient; Fieller et al., 1957
zmap = rho2z(diffmap);
zmap(zmap < 1.96) = 0;

% difference (subordinate + basic) vs. superordinate (i.e., [-0.5 -0.5 1])
diffmap2 = z2r((squeeze(r2z(rho(1,:,:)))*-0.5 + squeeze(r2z(rho(2,:,:)))*-0.5) + squeeze(r2z(rho(3,:,:)))*1);
zmap2 = rho2z(diffmap2);
zmap2(zmap2 < 1.96) = 0;

figure(3);
subplot(2,2,1)
imagesc([-0.6:0.02:1.18],30:-1:4, diffmap, [-0.06 0.06]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Spearman''s rho';
width=1200;
height=600;
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title('basic - subordinate');
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

subplot(2,2,2)
imagesc([-0.6:0.02:1.18],30:-1:4, zmap, [-3 3]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'z score';
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title('One-sided p = 0.025');
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

subplot(2,2,3)
imagesc([-0.6:0.02:1.18],30:-1:4, diffmap2, [-0.06 0.06]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Spearman''s rho';
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title('superordinate - (basic + subordinate)');
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

subplot(2,2,4)
imagesc([-0.6:0.02:1.18],30:-1:4, zmap2, [-3 3]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'z score';
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title('One-sided p = 0.025');
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

% randomization for cluster size correction
% H0: true difference between rhos is 0
% test cluster size for random permutations of RDM vector
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
mxc{nmr} = 0;
else
mxc{nmr} = [l.Area];
end
end
save([dapath, 'mxczmap.mat'], 'mxc');
threshed = zmap >= 1.96;
[cnum, N] = bwlabel(threshed, 4); % 4-connected (no vertical connection)
clusters    = arrayfun(@(x) find(cnum==x), 1:N, 'UniformOutput', false); % sorted labels 1:4
clusterSize = arrayfun(@(x) numel([x{:}]), clusters);


pSigClusters = arrayfun(@(x) sum([mxc{:}] > x) / numel([mxc{:}]), clusterSize);
sigIndex = find(pSigClusters < .05);
[rw,cl] = ind2sub([size(da_grandmean,3), size(da_grandmean,4)], clusters{2});
pltindex = [freqaxis(rw); timeaxis(cl)]'; % indices for plotting


