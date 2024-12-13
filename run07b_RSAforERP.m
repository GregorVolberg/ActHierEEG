% ==============
% get and plot RSA results
% ==============

% paths and files
addpath('./func');
dapath      = './DAerp/';
lowVisControlRDM = importdata('./Controlmodeldsm.mat'); % low level visual control RDM (72 x 72)
sceneControlRDM  = importdata('./places_last2nd.mat'); % Scene control RDM (72 x 72)
behavioralRDM    = importdata('./humanRating_mds.mat'); % behavioral RDM (12 x 12)

%% analysis 1: average DA in clusters, correlate with model RDMs
RDMs = get_modelRDMs();
da   = get_decodingaccuracy_fromfiles(dapath);

% average decoding accuracy per timebin per participant
nconditions = size(da.data,2);
ind = find(tril(ones(nconditions, nconditions), -1) == 1);
rho = nan(numel(RDMs), numel(da.subjects), size(da.data, 4));
for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
        for vp = 1:numel(da.subjects)
           for timebin = 1:size(da.data,4)
            tmpDA        =     squeeze(da.data(vp,:,:, timebin));
            rho(rdm, vp, timebin) = corr(RDM(ind), tmpDA(ind), 'type', 'Spearman'); % 3 x 22 x 600
            end
        end
end

% bootstrapping ci for rho 
r2z = @(x) atanh(x); % r to fisher's z
z2r = @(x) tanh(x);  % fisher's z to r
n_bootstrap = 1000;
rng(32); % initialize random number 
rnd               = arrayfun(@(x) randsample(size(rho, 2), size(rho, 2), 'true'), ...
                      1:n_bootstrap, 'UniformOutput', false);

%bs_rho  = arrayfun(@(x) z2r(shiftdim(squeeze(mean(r2z(rho(:,rnd{x},:)), 2))', 1)), ...
 %                      1:n_bootstrap, 'UniformOutput', false);
for x = 1:size(rnd,2)
    bs_rho(:,:,x) = z2r(squeeze(mean(r2z(rho(:,rnd{x},:)),2)));
end

[ci_lower(1,:), ci_upper(1,:)]  = deal(squeeze(quantile(bs_rho(1,:,:), 0.025, 3)), ...
                             squeeze(quantile(bs_rho(1,:,:), 0.975, 3))); 
[ci_lower(2,:), ci_upper(2,:)]  = deal(squeeze(quantile(bs_rho(2,:,:), 0.025, 3)), ...
                             squeeze(quantile(bs_rho(2,:,:), 0.975, 3))); 
[ci_lower(3,:), ci_upper(3,:)]  = deal(squeeze(quantile(bs_rho(3,:,:), 0.025, 3)), ...
                             squeeze(quantile(bs_rho(3,:,:), 0.975, 3))); 


% plot
cbPalette = {'#999999', '#E69F00', '#56B4E9', ...
               '009E73', '#F0E442', '#0072B2', ...
               '#D55E00', '#CC79A7'};
for j = 1:numel(cbPalette)
str = cbPalette{j};
cbcolor(j,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

figure(1);
plot(squeeze(mean(rho, 2))', 'LineWidth',1.5);
colororder(cbcolor);
hold on
fill([1:600, flip(1:600)], [ci_lower(3,:), flip(ci_upper(3,:))], 'r', 'FaceColor', cbPalette{3}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on
fill([1:600, flip(1:600)], [ci_lower(2,:), flip(ci_upper(2,:))], 'r', 'FaceColor', cbPalette{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on
fill([1:600, flip(1:600)], [ci_lower(1,:), flip(ci_upper(1,:))], 'r', 'FaceColor', cbPalette{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
set(gca,'XTick',[1:100:601], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([-0.05 0.15]);
xlabel('Time (s)'); ylabel('Correlation (Spearman''s rho)');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
rectFig = get(gcf,'position');
height=300;
width=height * 4/3;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
leg = legend('subordinate','basic', 'superordinate');
set(leg, 'box', 'off')


[mx, indx] = max(squeeze(mean(rho, 2))');
time = [-0.2:0.002:0.998];
conditions = {'subordinate','basic', 'superordinate'};
fprintf('\n\n');
for n = 1:numel(mx)
fprintf('%s:\t Maximum %.3f at %.3f s\n', conditions{n}, round(mx(n),3), time(indx(n)));
end

% % export data to ascii table for R
% nconditions = size(da.data,2);
% ind = find(tril(ones(nconditions, nconditions), -1) == 1);
% [r, c] = ind2sub([72,72], ind);
% testImg = zeros(72,72);
% testImg(ind) = ind;
% figure; imagesc(ind);
% images = strcat('R',num2str(r, '%02.f'), 'C', num2str(c, '%02.f'));

% standard RSA with betas
beta  = nan(numel(RDMs), numel(da.subjects), size(da.data, 4));
for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
        for vp = 1:numel(da.subjects)
           for timebin = 1:size(da.data,4)
            tmpDA        =     squeeze(da.data(vp,:,:, timebin));
            beta(rdm, vp, timebin) = regress(zscore(tmpDA(ind)), zscore(RDM(ind)));
            end
        end
end

figure(2);
plot(squeeze(z2r(mean(r2z(beta), 2)))', 'LineWidth',1.5);
colororder(cbcolor);

% same with average DA
betaAVG  = nan(numel(RDMs), size(da.data, 4));
daAVG    = squeeze(mean(da.data,1));
for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
           for timebin = 1:size(daAVG, 3)
            tmpDA        =     squeeze(daAVG(:,:, timebin));
            model        = fitlm(zscore(tmpDA(ind)), zscore(RDM(ind)));
            betaAVG(rdm, timebin) = model.Coefficients{'x1', 'Estimate'};
            
            ci = coefCI(model, 0.05);
            
            end
end

figure(3);
plot(betaAVG', 'LineWidth',1.5);
colororder(cbcolor);
k.coefCI

% regressing out the effect of scene and low-level-vision
% i.e. remove variance of covariates from criterion and predictor, correlate residuals
X        = [lowVisControlRDM(ind), ...
           sceneControlRDM(ind)];
betasc  = nan(numel(RDMs), numel(da.subjects), size(da.data, 4));
for rdm = 1:numel(RDMs)
    RDM = RDMs{rdm};
        for vp = 1:numel(da.subjects)
           for timebin = 1:size(da.data,4)
            tmpDA        =     squeeze(da.data(vp,:,:, timebin));
            model1       = fitlm(X, tmpDA(ind));
            model2       = fitlm(X, RDM(ind));
            betasc(rdm, vp, timebin) = regress(zscore(model1.Residuals.Raw), zscore(model2.Residuals.Raw));
            end
        end
end

figure(4);
plot(squeeze(z2r(mean(r2z(betasc), 2)))', 'LineWidth',1.5);
colororder(cbcolor);




