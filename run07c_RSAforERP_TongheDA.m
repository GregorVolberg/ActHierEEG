% ==============
% get and plot RSA results
% ==============
% use Tonghes decoding accuracy data and model RDMs
% paths and files
addpath('./func');
%dapath      = './DAerp/';
lowVisControlRDM = importdata('./Controlmodeldsm.mat'); % low level visual control RDM (72 x 72)
sceneControlRDM  = importdata('./places_last2nd.mat'); % scene control RDM (72 x 72)
behavioralRDM    = importdata('./humanRating_mds.mat'); % behavioral RDM (12 x 12)

%% analysis 1: per-participant DA, correlate with model RDMs
% get RDMS and decoding accoracy
RDMs = {};
RDMs{1} = importdata('./osf_travelfmri/sub.mat');
RDMs{2} = importdata('./osf_travelfmri/basic.mat');
RDMs{3} = importdata('./osf_travelfmri/super.mat');

tzpath = '../Results/RSA/EEG_based_RDM/';
tmp = dir([tzpath, '*_decAll.mat']);
tmpchar = char({tmp.name});
da.files = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
da.subjects  = cellstr(tmpchar(:,1:5));
da.data = deal(nan(22, 72, 72, 600));
for vp = 1:numel(da.subjects)
    tmp = load(da.files{vp});
    da.data(vp, :,:,:) = tmp.res.diss; 
end

% get indices for lower triangle of RDMs and decoding accuracy matrices
%nconditions = size(da.data,2);
%ind = find(tril(ones(nconditions, nconditions), -1) == 1);

% get predictors for multiple regression analysis
X        = [squareform(lowVisControlRDM); ...
           squareform(sceneControlRDM)]';

% pre-allocate results matrices (3 x 22 x 600) 
[rho, beta, betam, betat] = deal(nan(numel(RDMs), numel(da.subjects), size(da.data, 4)));

% loop over RDMs, participants and time bins
% use squareform instead of indicator
for rdm = 1:numel(RDMs)
    RDM = squareform(RDMs{rdm})';
        for vp = 1:numel(da.subjects)
           for timebin = 1:size(da.data,4)
            tmpDA  =     squareform(squeeze(da.data(vp,:,:, timebin)))';
            tongheDA = squareform(EEG_reshape_to_fMRI(squeeze(da.data(vp,:,:, timebin))))'; % tonghe model
            rho(rdm, vp, timebin) = corr(RDM, tmpDA, 'type', 'Spearman');
            model0 = fitlm(zscore(RDM), zscore(tmpDA));
            modelt = fitlm(zscore(RDM), zscore(tongheDA)); % tonghe model
            beta(rdm, vp, timebin) = model0.Coefficients{'x1', 'Estimate'};
            betat(rdm, vp, timebin) = modelt.Coefficients{'x1', 'Estimate'};
            model1 = fitlm(X, tmpDA);
            model2 = fitlm(X, RDM);
            model3 = fitlm(zscore(model1.Residuals.Raw), zscore(model2.Residuals.Raw));
            betam(rdm, vp, timebin)  = model3.Coefficients{'x1', 'Estimate'};
            end
        end
end

% get bootstrapping CI's
r2z = @(x) atanh(x); % r to fisher's z
z2r = @(x) tanh(x);  % fisher's z to r
n_bootstrap = 1000;
rng(32); % initialize random number 
rnd  = arrayfun(@(x) randsample(size(rho, 2), size(rho, 2), 'true'), ...
               1:n_bootstrap, 'UniformOutput', false);

rhoCI   = get_ci(rho, rnd);
betaCI  = get_ci(beta, rnd);
betamCI = get_ci(betam, rnd);
betatCI = get_ci(betat, rnd);

%% plotting
cbcolor = get_cbPalette;

% Spearman
figure(1);
plot(squeeze(z2r(mean(r2z(betat), 2)))', 'LineWidth',1.5);
colororder(cbcolor);
hold on; fill([1:600, flip(1:600)], [squeeze(betatCI(:, 3, 1))', flip(squeeze(betatCI(:, 3, 2)))'], 'r', 'FaceColor', cbcolor(3, :), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(betatCI(:, 2, 1))', flip(squeeze(betatCI(:, 2, 2)))'], 'r', 'FaceColor', cbcolor(2, :), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(betatCI(:, 1, 1))', flip(squeeze(betatCI(:, 1, 2)))'], 'r', 'FaceColor', cbcolor(1, :), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
sigtime = (rhoCI(:,:,1) > 0) .* [-0.02 -0.03 -0.04];
sigtime(sigtime == 0) = nan;
hold on; plot(sigtime(:,1), 'Linewidth', 1.5, 'Color', cbcolor(1,:));
hold on; plot(sigtime(:,2), 'Linewidth', 1.5, 'Color', cbcolor(2,:));
hold on; plot(sigtime(:,3), 'Linewidth', 1.5, 'Color', cbcolor(3,:));
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.05 0.15]);
xlabel('Time (s)'); ylabel('beta');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')

% beta
figure(2);
plot(squeeze(z2r(mean(r2z(beta), 2)))', 'LineWidth',1.5);
colororder(cbcolor);
hold on; fill([1:600, flip(1:600)], [squeeze(betaCI(:, 3, 1))', flip(squeeze(betaCI(:, 3, 2)))'], 'r', 'FaceColor', cbcolor(3,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(betaCI(:, 2, 1))', flip(squeeze(betaCI(:, 2, 2)))'], 'r', 'FaceColor', cbcolor(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(betaCI(:, 1, 1))', flip(squeeze(betaCI(:, 1, 2)))'], 'r', 'FaceColor', cbcolor(1,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
sigtime = (betaCI(:,:,1) > 0) .* [-0.02 -0.03 -0.04];
sigtime(sigtime == 0) = nan;
hold on; plot(sigtime(:,1), 'Linewidth', 1.5, 'Color', cbcolor(1,:));
hold on; plot(sigtime(:,2), 'Linewidth', 1.5, 'Color', cbcolor(2,:));
hold on; plot(sigtime(:,3), 'Linewidth', 1.5, 'Color', cbcolor(3,:));
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.05 0.15]);
xlabel('Time (s)'); ylabel('beta');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')

% betam
figure(3);
plot(squeeze(z2r(mean(r2z(betam), 2)))', 'LineWidth',1.5);
colororder(cbcolor);
ci = betamCI;
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 3, 1))', flip(squeeze(ci(:, 3, 2)))'], 'r', 'FaceColor', cbcolor(3,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 2, 1))', flip(squeeze(ci(:, 2, 2)))'], 'r', 'FaceColor', cbcolor(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 1, 1))', flip(squeeze(ci(:, 1, 2)))'], 'r', 'FaceColor', cbcolor(1,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
sigtime = (ci(:,:,1) > 0) .* [-0.02 -0.03 -0.04];
sigtime(sigtime == 0) = nan;
hold on; plot(sigtime(:,1), 'Linewidth', 1.5, 'Color', cbcolor(1,:));
hold on; plot(sigtime(:,2), 'Linewidth', 1.5, 'Color', cbcolor(2,:));
hold on; plot(sigtime(:,3), 'Linewidth', 1.5, 'Color', cbcolor(3,:));
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.05 0.15]);
xlabel('Time (s)'); ylabel('beta control model');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')

%% same with average DA
% pre-allocate results matrices (3 x 22 x 600) 
[rhoAVG, betaAVG, betamAVG] = deal(nan(numel(RDMs), size(da.data, 4)));
daAVG    = squeeze(mean(da.data,1));
clear tmpDA
% loop over RDMs and time bins
for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
       for timebin = 1:size(daAVG,3)
            tmpDA  = squeeze(daAVG(:,:, timebin));
            rhoAVG(rdm, timebin) = corr(RDM(ind), tmpDA(ind), 'type', 'Spearman');
            model0 = fitlm(zscore(RDM(ind)), zscore(tmpDA(ind)));
            betaAVG(rdm, timebin) = model0.Coefficients{'x1', 'Estimate'};
            tmp = coefCI(model0, 0.05);
            ciAVG(timebin,rdm, :) = tmp(2,:);
            model1 = fitlm(X, tmpDA(ind));
            model2 = fitlm(X, RDM(ind));
            model3 = fitlm(zscore(model1.Residuals.Raw), zscore(model2.Residuals.Raw));
            betamAVG(rdm, timebin)  = model3.Coefficients{'x1', 'Estimate'};
            tmp = coefCI(model3, 0.05);
            cimAVG(timebin,rdm, :) = tmp(2,:);
       end
end

%% plot
% AVG rho
figure(4);
plot(rhoAVG', 'LineWidth',1.5);
colororder(cbcolor);
% ci = rhoCI;
% hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 3, 1))', flip(squeeze(ci(:, 3, 2)))'], 'r', 'FaceColor', cbPalette{3}, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.2);
% hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 2, 1))', flip(squeeze(ci(:, 2, 2)))'], 'r', 'FaceColor', cbPalette{2}, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.2);
% hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 1, 1))', flip(squeeze(ci(:, 1, 2)))'], 'r', 'FaceColor', cbPalette{1}, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.2);
% sigtime = (ci(:,:,1) > 0) .* [-0.02 -0.03 -0.04];
% sigtime(sigtime == 0) = nan;
% hold on; plot(sigtime(:,1), 'Linewidth', 1.5, 'Color', cbcolor(1,:));
% hold on; plot(sigtime(:,2), 'Linewidth', 1.5, 'Color', cbcolor(2,:));
% hold on; plot(sigtime(:,3), 'Linewidth', 1.5, 'Color', cbcolor(3,:));
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.1 0.35]);
xlabel('Time (s)'); ylabel('Average rho');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')

% avg beta
figure(5);
plot(betaAVG', 'LineWidth',1.5);
colororder(cbcolor);
ci = ciAVG;
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 3, 1))', flip(squeeze(ci(:, 3, 2)))'], 'r', 'FaceColor', cbcolor(3,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 2, 1))', flip(squeeze(ci(:, 2, 2)))'], 'r', 'FaceColor', cbcolor(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 1, 1))', flip(squeeze(ci(:, 1, 2)))'], 'r', 'FaceColor', cbcolor(1,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
sigtime = (ci(:,:,1) > 0) .* [-0.02 -0.03 -0.04];
sigtime(sigtime == 0) = nan;
hold on; plot(sigtime(:,1), 'Linewidth', 1.5, 'Color', cbcolor(1,:));
hold on; plot(sigtime(:,2), 'Linewidth', 1.5, 'Color', cbcolor(2,:));
hold on; plot(sigtime(:,3), 'Linewidth', 1.5, 'Color', cbcolor(3,:));
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.1 0.35]);
xlabel('Time (s)'); ylabel('Average beta');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')

% avg beta control
figure(6);
plot(betamAVG', 'LineWidth',1.5);
colororder(cbcolor);
ci = cimAVG;
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 3, 1))', flip(squeeze(ci(:, 3, 2)))'], 'r', 'FaceColor', cbcolor(3,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 2, 1))', flip(squeeze(ci(:, 2, 2)))'], 'r', 'FaceColor', cbcolor(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on; fill([1:600, flip(1:600)], [squeeze(ci(:, 1, 1))', flip(squeeze(ci(:, 1, 2)))'], 'r', 'FaceColor', cbcolor(1,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
sigtime = (ci(:,:,1) > 0) .* [-0.02 -0.03 -0.04];
sigtime(sigtime == 0) = nan;
hold on; plot(sigtime(:,1), 'Linewidth', 1.5, 'Color', cbcolor(1,:));
hold on; plot(sigtime(:,2), 'Linewidth', 1.5, 'Color', cbcolor(2,:));
hold on; plot(sigtime(:,3), 'Linewidth', 1.5, 'Color', cbcolor(3,:));
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.1 0.35]);
xlabel('Time (s)'); ylabel('Average beta control model');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')
% 
% 
% [mx, indx] = max(squeeze(mean(rho, 2))');
% time = [-0.2:0.002:0.998];
% conditions = {'subordinate','basic', 'superordinate'};
% fprintf('\n\n');
% for n = 1:numel(mx)
% fprintf('%s:\t Maximum %.3f at %.3f s\n', conditions{n}, round(mx(n),3), time(indx(n)));
% end
% 
% 
% 
% 
% % regressing out the effect of scene and low-level-vision
% % i.e. remove variance of covariates from criterion and predictor, correlate residuals
% X        = [lowVisControlRDM(ind), ...
%            sceneControlRDM(ind)];
% betasc  = nan(numel(RDMs), numel(da.subjects), size(da.data, 4));
% for rdm = 1:numel(RDMs)
%     RDM = RDMs{rdm};
%         for vp = 1:numel(da.subjects)
%            for timebin = 1:size(da.data,4)
%             tmpDA        =     squeeze(da.data(vp,:,:, timebin));
%             model1       = fitlm(X, tmpDA(ind));
%             model2       = fitlm(X, RDM(ind));
%             betasc(rdm, vp, timebin) = regress(zscore(model1.Residuals.Raw), zscore(model2.Residuals.Raw));
%             end
%         end
% end
% 
% figure(4);
% plot(squeeze(z2r(mean(r2z(betasc), 2)))', 'LineWidth',1.5);
% colororder(cbcolor);
% 
% 
% 
% 
