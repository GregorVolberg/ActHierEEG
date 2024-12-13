% ==============
% get and plot RSA results
% ==============

%% paths and files
addpath('./func');
dapath      = './DA/';
clusterfile = [dapath, 'clustermat.mat'];
rhobootstrapfile = [dapath, 'rho_per_bin.mat'];  

%% analysis 1: average DA in clusters, correlate with model RDMs
RDMs = get_modelRDMs();
da   = get_decodingaccuracy_fromfiles(dapath);
load(clusterfile); % contains clustermat

% average decoding accuracy per cluster per participant
da_per_cluster = arrayfun(@(x) mean(mean(da.data(:,:,:,clustermat(x).cluster(:,1:2)),5),4), ...
                 1:2, 'UniformOutput', false); % 2 cells with 22 x 72 x 72 matrices

% correlate lower triangle data of decoding accuracy and RDMs 
nconditions = size(da.data,2);
ind = find(tril(ones(nconditions, nconditions), -1) == 1);
rho = nan(numel(clustermat), numel(RDMs), numel(da.subjects));
for clster = 1:numel(clustermat)
    da_cluster = da_per_cluster{clster};
    for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
        for vp = 1:numel(da.subjects)
        daclusterPP =     squeeze(da_cluster(vp,:,:));
        rho(clster, rdm, vp) = corr (RDM(ind), daclusterPP(ind), 'type', 'Spearman'); % 2 x 3 x 22
        end
    end
end
             
% 
% nconditions = size(da.data,2);
% c1mean = nan(22, nconditions, nconditions);
% c2mean = nan(22, nconditions, nconditions);
% for vp = 1:numel(dafiles.vps)
%     load(dafiles.name{vp});
%     for con1 = 1:size(DAmean, 1)
%         for con2 = 1:size(DAmean, 2)
%         tmpmat = squeeze(DAmean(con1, con2, :, :));
%         c1mean(vp, con1, con2) = mean(tmpmat(clustermat(1).cluster(:,3)));
%         c2mean(vp, con1, con2) = mean(tmpmat(clustermat(2).cluster(:,3)));
%         end
%     end
% end
% 

% rho1m = mean(rho1,2); rho1sd = std(rho1'); rho1se = rho1sd / sqrt(size(rho1,2));
% rho2m = mean(rho2,2); rho2sd = std(rho2'); rho2se = rho2sd / sqrt(size(rho2,2));
% ycoord = [rho1m-rho1se', rho1m+rho1se'; rho2m-rho2se', rho2m+rho2se']';
% xcoord = [repmat(1:3,2,1), repmat(5:7,2,1)];

rho1m = mean(rho(1,:,:),3); rho1se = std(squeeze(rho(1,:,:))') / sqrt(numel(da.subjects));
rho2m = mean(rho(2,:,:),3); rho2se = std(squeeze(rho(2,:,:))') / sqrt(numel(da.subjects));
ycoord = [rho1m-rho1se', rho1m+rho1se'; rho2m-rho2se', rho2m+rho2se']';
xcoord = [repmat(1:3,2,1), repmat(5:7,2,1)];

% noise ceiling
% To estimate the upper bound we correlated (Spearman’s R) each participant’s neural RDM
% with the mean neural RDM across all participants. 
% To estimate the lower bound we correlated (Spearman’s R) each participant’s neural RDM 
% with the mean neural RDM excluding that participant iteratively for all participants. 
% We averaged the results, yielding estimates of the lower and upper noise ceiling for infants and adults
% upperbound
rho1ceil = nan(size(c1mean,1), 2);
rho2ceil = nan(size(c1mean,1), 2);

for vp = 1:size(c1mean, 1)
    c1      = squeeze(c1mean(vp, :,:));
    c1upper = squeeze(mean(c1mean, 1));
    c1l     = c1mean; 
    c1l(vp,:,:) = [];
    c1lower = squeeze(mean(c1l, 1));
    
    c2      = squeeze(c2mean(vp, :,:));
    c2upper = squeeze(mean(c2mean, 1));
    c2l     = c2mean; 
    c2l(vp,:,:) = [];
    c2lower = squeeze(mean(c2l, 1));
    
    rho1ceil(vp, 1) = corr (c1(ind), c1lower(ind), 'type', 'Spearman');
    rho1ceil(vp, 2) = corr (c1(ind), c1upper(ind), 'type', 'Spearman');
    rho2ceil(vp, 1) = corr (c2(ind), c2lower(ind), 'type', 'Spearman');
    rho2ceil(vp, 2) = corr (c2(ind), c2upper(ind), 'type', 'Spearman');
end

noiseceiling1 = mean(rho1ceil);
noiseceiling2 = mean(rho2ceil);

% plot
figure(1);
jitterwidth = 0.1; % uniform -0.1:0.1
lightblue = [0.6784    0.8471    0.9020];
x = reshape(repmat([1, 2, 3, 5, 6, 7], size(rho1, 2),1), [size(rho1, 2)*6, 1]);
x = add_jitter(x, jitterwidth);
y = [reshape(rho1', [prod(size(rho1)), 1]); reshape(rho2', [prod(size(rho2)), 1])];
bar([1, 2, 3, 5, 6, 7], [rho1m, rho2m], 0.4, ...
    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
rectangle('Position', [0.3 noiseceiling1(1) 3.4 noiseceiling1(2)-noiseceiling1(1)], ...
    'FaceColor',[0.9 0.9 0.9], 'LineWidth', 0.1, 'EdgeColor', [0.9 0.9 0.9]);
rectangle('Position', [4.3 noiseceiling2(1) 3.4 noiseceiling2(2)-noiseceiling2(1)], ...
    'FaceColor',[0.9 0.9 0.9], 'LineWidth', 0.1, 'EdgeColor', [0.9 0.9 0.9]);
hold on; scatter(x, y, 10, lightblue, 'filled');
hold on; line(xcoord, ycoord, 'Linewidth', 2, 'Color', 'k');
ylim([-0.04 0.30]);
xlim([0.2 7.7]);
set(gca,'XTick',[1:3, 5:7], ...
    'XTickLabel', {'subordinate', 'basic', 'superordinate'}, ...
    'YTick',[0:0.04:0.28], ...
    'YTickLabel', [0:0.04:0.28]);
xlabel('Level'); ylabel('Correlation (Spearman''s rho)');
text(2, -0.02, 'delta cluster', 'HorizontalAlignment','center');
text(6, -0.02, 'alpha cluster', 'HorizontalAlignment','center');
title('Model RDMs versus decoding accuracy');
rectFig = get(gcf,'position');
width=600;
height=300;
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

% bootstrapping ci 

for btstrp = 1:1000
rnd = arrayfun(@(x) randsample(size(rho1, 2), size(rho1, 2), 'true'), 1:3, 'UniformOutput', false);
mm  = arrayfun(@(x) mean(rho1(x, rnd{x})), 1:3);
bs(btstrp,:) = mm;
end
quantile(bs, 0.025)
quantile(bs, 0.975) % all significant

for btstrp = 1:1000
rnd = arrayfun(@(x) randsample(size(rho2, 2), size(rho2, 2), 'true'), 1:3, 'UniformOutput', false);
mm  = arrayfun(@(x) mean(rho2(x, rnd{x})), 1:3);
bs2(btstrp,:) = mm;
end
quantile(bs2, 0.025)
quantile(bs2, 0.975) % all significant


% per TF bin correlation
% mean vp
da_grandmean = squeeze(mean(da.data, 1));
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
CH = colorbar('eastoutside'); CH.Label.String = 'Correlation (Spearman)';
title(titles{j});
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
end
rectFig = get(gcf,'position');
width=600;
height=900;
set(gcf,'position',[rectFig(1),rectFig(2)- height,width,height], 'color', 'white');

%find max
col={'k', 'b', 'r'};
timeaxis = [-0.6:0.02:1.18];
freqaxis = [4:30];

figure(2);
[m]      = arrayfun(@(x)  max(squeeze(rho(x,:,:))), 1:3, 'UniformOutput', false); % time
[mt,it]  = arrayfun(@(x) max(x{:}), m)
[m2]     = arrayfun(@(x)  max(squeeze(rho(x,:,:))'), 1:3, 'UniformOutput', false); % freq
[m2t,i2t]= arrayfun(@(x) max(x{:}), m2);

plot(timeaxis, cell2mat(m')', 'LineWidth', 1.5);
title('Maximum Correlation RSA');
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
xlabel('Time (ms)'); ylabel('Correlation (Spearman)');
rectFig = get(gcf,'position');
width=600;
height=300;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
legend({'subordinate', 'basic', 'superordinate'})
legend('boxoff')
arrayfun(@(x, y) text(timeaxis(x), y+0.015, num2str(timeaxis(x)), 'HorizontalAlignment', 'center'), it, mt)

plot(cell2mat(m2')'*-1, freqaxis, 'LineWidth', 1.5);
title('Maximum Correlation RSA');
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
xlabel('Correlation (Spearman)'); ylabel('Frequency (Hz)');
ylim([4 30]); xlim([-0.2 0.2])
set(gca,'YTick',[5:5:30], ...
    'YTickLabel', [5:5:30]);
set(gca,'XTick',[-0.2:0.05:0.2], ...
    'XTickLabel', [0.2:-0.05:-0.2]);
rectFig = get(gcf,'position');
width=300;
height=600;
set(gcf,'position',[rectFig(1),rectFig(2)-height/2,width,height], 'color', 'white');
legend({'subordinate', 'basic', 'superordinate'})
legend('boxoff')
arrayfun(@(x, y) text(-x+0.02, freqaxis(y), num2str(freqaxis(y)), 'VerticalAlignment', 'middle'), m2t, i2t)

% CI
load(rhobootstrapfile); % contains struct rho_per_bin
log=rho_per_bin.rho_boot > shiftdim(repmat(rho_per_bin.rho, [1,1,1,1000]),3);
log2=squeeze(sum(log,1))./size(rho_per_bin.rho_boot,1);
imagesc(squeeze(1-log2(1,:,:)), [0.9 1])



