% ==============
% get and plot RSA results, per cluster
% ==============

%% paths and files
addpath('./func');
dapath      = './DA/';
clusterfile = [dapath, 'clustermat.mat'];

%% analysis 1: average DA in clusters, correlate with model RDMs
RDMs = get_modelRDMs();
da   = get_decodingaccuracy_fromfiles(dapath);
load(clusterfile); % contains clustermat

% average decoding accuracy per cluster per participant
nconditions = size(da.data,2);
cmean = nan(2, 22, nconditions, nconditions); % cluster mean
for vp = 1:numel(da.files)
    for con1 = 1:size(da.data, 2)
        for con2 = 1:size(da.data, 3)
        tmpmat = squeeze(da.data(vp, con1, con2, :, :));
        cmean(1, vp, con1, con2) = mean(tmpmat(clustermat(1).cluster(:,3)));
        cmean(2, vp, con1, con2) = mean(tmpmat(clustermat(2).cluster(:,3)));
        end
    end
end

ind = find(tril(ones(nconditions, nconditions), -1) == 1);
rho = nan(numel(clustermat), numel(RDMs), numel(da.subjects));
for clster = 1:numel(clustermat)
    da_cluster = squeeze(cmean(clster,:,:,:));
    for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
        for vp = 1:numel(da.subjects)
        daclusterPP =     squeeze(da_cluster(vp,:,:));
        rho(clster, rdm, vp) = corr (RDM(ind), daclusterPP(ind), 'type', 'Spearman'); % 2 x 3 x 22
        end
    end
end

% bootstrapping ci for rho 
r2z = @(x) atanh(x); % r to fisher's z
z2r = @(x) tanh(x);  % fisher's z to r
n_bootstrap = 1000;
rng(32); % initialize random number 
rnd               = arrayfun(@(x) randsample(size(rho, 3), size(rho, 3), 'true'), ...
                      1:n_bootstrap, 'UniformOutput', false);
bootstrapped_rho  = arrayfun(@(x) z2r(reshape(squeeze(mean(r2z(rho(:,:, rnd{x})),3))', 1,6)), ...
                      1:n_bootstrap, 'UniformOutput', false);
bs_rho = cell2mat(bootstrapped_rho');
sig    = find(quantile(bs_rho, 0.025) > 0); % CIs do not include 0, -> all significant

% ANOVA for cluster x model
aovdat = cell2mat(arrayfun(@(x, y) squeeze(r2z(rho(x,y,:))), ...
            [1 1 1 2 2 2], [1 2 3 1 2 3], 'UniformOutput', false));
data = array2table([[1:size(aovdat,1)]', aovdat], 'VariableNames', ...
         {'id', 'c1_r1', 'c1_r2', 'c1_r3',  'c2_r1', 'c2_r2', 'c2_r3'});
w = table(categorical([1 1 1 2 2 2].'), categorical([1 2 3 1 2 3].'), ...
         'VariableNames', {'cluster', 'rdm'}); % within-desing
rm  = fitrm(data, 'c1_r1-c2_r3 ~ 1', 'WithinDesign', w);
res = ranova(rm, 'withinmodel', 'cluster*rdm');
disp(res);

% post-hoc I: main effect RDM
fprintf('\n\n'); 
mainRDM = (aovdat(:, 1:3) + aovdat(:, 4:6))/2;
[~, p, ~, t] = ttest(mainRDM(:,1) - mainRDM(:,2)); % subordinate vs. basic
fprintf('subordinate vs. basic, t(%i) = %.3f, p = %.3f\n', t.df, t.tstat, p);
[~, p, ~, t] = ttest(mainRDM(:,1) - mainRDM(:,3)); % subordinate vs. superordinate
fprintf('subordinate vs. superordinate, t(%i) = %.3f, p = %.3f\n', t.df, t.tstat, p);
[~, p, ~, t] = ttest(mainRDM(:,2) - mainRDM(:,3)); % basic vs. superordinate
fprintf('basic vs. superordinate, t(%i) = %.3f, p = %.3f\n', t.df, t.tstat, p);

% post-hoc II: interaction cluster x rdm
conds = [1 4; 2 5; 3 6]; 
condnames = {'subordinate delta', 'basic delta', 'superordinate delta', ...
             'subordinate alpha', 'basic alpha', 'superordinate alpha'};
condmeans = mean(aovdat);
fprintf('\n\n'); 
for j = 1:3
[~,p,~,t] = ttest(aovdat(:,conds(j,1)) - aovdat(:, conds(j,2)));
fprintf('%s  vs. %s, %.3f vs. %.3f, t(%i) = %.3f, p = %.3f\n', ...
    condnames{conds(j,1)}, condnames{conds(j,2)}, condmeans(conds(j,1)), condmeans(conds(j,2)), ...
    t.df, t.tstat, p);
end
fprintf('\n\n'); 

% coordinates for plotting
rho1m = z2r(mean(r2z(rho(1,:,:)),3)); rho1se = z2r(std(squeeze(r2z(rho(1,:,:)))') / sqrt(numel(da.subjects)));
rho2m = z2r(mean(r2z(rho(2,:,:)),3)); rho2se = z2r(std(squeeze(r2z(rho(2,:,:)))') / sqrt(numel(da.subjects)));
ycoord = [rho1m-rho1se, rho2m-rho2se; rho1m+rho1se, rho2m+rho2se];
xcoord = [repmat(1:3,2,1), repmat(5:7,2,1)];
xcoordsig = xcoord(1, sig);
ycoordsig = ycoord(2,sig) + 0.02; % place slightly above SE bar

% noise ceiling
% To estimate the upper bound we correlated (Spearman’s R) each participant’s neural RDM
% with the mean neural RDM across all participants. 
% To estimate the lower bound we correlated (Spearman’s R) each participant’s neural RDM 
% with the mean neural RDM excluding that participant iteratively for all participants. 
% We averaged the results, yielding estimates of the lower and upper noise ceiling for infants and adults
% upperbound
noiseceil = nan(2, size(cmean,1), 2); % cluster x subject x [lower, upper]
for clstr = 1:2
for vp = 1:size(cmean, 2)
    vpdat    = squeeze(cmean(clstr,vp, :,:));
    nc_upper = squeeze(mean(cmean(clstr,:,:,:), 2));
    tmpmat   = squeeze(cmean(clstr,:,:,:)); tmpmat(vp,:,:) = [];
    nc_lower = squeeze(mean(tmpmat, 1));
    noiseceil(clstr, vp, 1) = corr (vpdat(ind), nc_lower(ind), 'type', 'Spearman');
    noiseceil(clstr, vp, 2) = corr (vpdat(ind), nc_upper(ind), 'type', 'Spearman');
end
end
noiseceiling = squeeze(mean(noiseceil,2));

% plot
figure(1);
jitterwidth = 0.1; % uniform -0.1:0.1
lightblue = [0.6784    0.8471    0.9020];
x = reshape(repmat([1, 2, 3, 5, 6, 7], size(rho, 3),1), [size(rho, 3)*6, 1]);
x = add_jitter(x, jitterwidth);
rhoplot = shiftdim(rho,2);
y = [reshape(squeeze([rhoplot(:,1,1:3)]), 22*3, 1); ...
    reshape(squeeze([rhoplot(:,2,1:3)]), 22*3, 1)]; 
bar([1, 2, 3, 5, 6, 7], [rho1m, rho2m], 0.4, ...
    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
rectangle('Position', [0.3 noiseceiling(1,1) 3.4 noiseceiling(1,2)-noiseceiling(1,1)], ...
    'FaceColor',[0.9 0.9 0.9], 'LineWidth', 0.1, 'EdgeColor', [0.9 0.9 0.9]);
rectangle('Position', [4.3 noiseceiling(2,1) 3.4 noiseceiling(2,2)-noiseceiling(2,1)], ...
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
text(xcoordsig, ycoordsig, '*', 'HorizontalAlignment','center', 'FontSize',14);
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
%stats
%line([1, 5], [0.13 0.13], 'LineStyle', '-', 'LineWidth', 1, 'Color', 'k');
line([2, 6], [0.16 0.16], 'LineStyle', '-', 'LineWidth', 1, 'Color', 'k');
%line([3, 7], [0.19 0.19], 'LineStyle', '-', 'LineWidth', 1, 'Color', 'k');
text(4, 0.18, 'p = .008', 'HorizontalAlignment','center', 'FontSize',12);
