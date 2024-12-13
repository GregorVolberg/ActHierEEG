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

% find maximum rho per rdm
timeaxis = [-0.6:0.02:1.18];
freqaxis = [4:30];
freqaxisplot = [30:-1:4];
[~, itime] = arrayfun(@(x) max(max(squeeze(rho(x,:,:)))), 1:3); % time points
[m, ifreq] = arrayfun(@(x) max(max(squeeze(rho(x,:,:))')), 1:3); % frequencies

% make TF plots with Spearman's rho per bin
titles = {'Subordinate Level', 'Basic Level', 'Superordinate Level'};
scales = {[0,0.15], [0,0.15], [0,0.15]};
figure(1);
for j = 1:3
fprintf('Spearmans''s rho = %.3f at %i Hz and %.3f\n', m(j), freqaxis(ifreq(j)), timeaxis(itime(j)));
subplot(3,1,j);
imagesc([-0.6:0.02:1.18],30:-1:4, squeeze(rho(j,:,:)), scales{j});
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Spearman''s rho';
title(titles{j});
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
text(timeaxis(itime(j)), freqaxisplot(ifreq(j)), '+', 'FontSize', 14, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center', 'Color', 'r')
end
rectFig = get(gcf,'position');
width=600;
height=900;
set(gcf,'position',[rectFig(1),rectFig(2)- height,width,height], 'color', 'white');

%% plot rho and CI over frequencies at time point of maximum 
alph  = 0.01;
zcrit = norminv(1-(alph/2)); % find ci from z score 
n     = length(ind);
spearmanCI = @(x) tanh(atanh(x) + zcrit * [-1, 1] / sqrt(1.06*(n-3))); % https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation incl. comments
figure(2);
for rdm = 1:3
rs = rho(rdm, :, itime(rdm));
ci = spearmanCI(rs');
subplot(3,1,rdm);
poly = polyshape([[-1*ci(:,1)'; 4:30]'; [-1*flipud(ci(:,2))'; 30:-1:4]']);
plot(poly, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]); hold on
line([0 0], [4 30], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.7 0.7 0.7]);
plot(rs*-1, freqaxis, 'LineWidth', 1.5, 'Color', 'b');
xlabel('Spearman''s rho'); ylabel('Frequency (Hz)');
ylim([4 30]); xlim([-0.21 0.21])
set(gca,'YTick',[5:5:30], ...
    'YTickLabel', [5:5:30]);
set(gca,'XTick',[-0.2:0.1:0.2], ...
    'XTickLabel', {'0.2', '0.1', '0.0', '-0.1', '-0.2'});%
width=200;
height=900;
line([0 -m(rdm)], [freqaxis(ifreq(rdm)), freqaxis(ifreq(rdm))], 'Color', 'r'); % ,maximum
text(-rs(ci(:,1) > 0), freqaxis(ci(:,1) > 0)-0.5, '*', 'FontSize', 18, 'VerticalAlignment', 'middle',...
                            'HorizontalAlignment', 'center', 'Color', 'b');
end
rectFig = get(gcf,'position');
set(gcf,'position',[rectFig(1),rectFig(2)-height/2,width,height], 'color', 'white');
 
%% plot rho and CI over time points at frequency of maximum 
figure(3);   
for rdm = 1:3
rs = squeeze(rho(rdm, ifreq(rdm), :));
ci = spearmanCI(rs);
subplot(3,1,rdm);
poly = polyshape([[timeaxis; ci(:,1)']'; [fliplr(timeaxis); flipud(ci(:,2))']']);
plot(poly, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]); hold on
line([0 0], [-0.1 0.2], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.7 0.7 0.7]);
line([-0.6 1.2], [0 0], 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.7 0.7 0.7]);
plot(timeaxis, rs, 'LineWidth', 1.5, 'Color', 'b');
title([num2str(freqaxis(ifreq(rdm))), ' Hz']);
xlabel('Time (s)'); ylabel('Spearman''s rho');
ylim([-0.1 0.21]); xlim([-0.6 1.2])
set(gca,'XTick',[-0.6:0.2:1.0], ...
    'XTickLabel', {'-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6', '0.8', '1.0'});
line([timeaxis(itime(rdm)), timeaxis(itime(rdm))], [0 m(rdm)], 'Color', 'r'); 
text(timeaxis(ci(:,1) > 0), rs(ci(:,1) > 0)-0.005, '*', 'FontSize', 18, 'VerticalAlignment', 'middle',...
                            'HorizontalAlignment', 'center', 'Color', 'b');
end
rectFig = get(gcf,'position');
width=600;
height=900;
set(gcf,'position',[rectFig(1),rectFig(2)-height/2,width,height], 'color', 'white');

%% Difference maps and cluster correction
rho_conditions  = [2 1; ...  
                   2 3; ...  
                   1 3];     
rho_labels = {'basic minus subordinate';  ...
              'basic minus superordinate'; ... 
              'subordinate minus superordinate'};
               
pn = [1, -1];
load([dapath, 'mxczmap.mat']);

figure(4);

for rdms = 1:3
diffmap = z2r(squeeze(r2z(rho(rho_conditions(rdms, 1),:,:))) - squeeze(r2z(rho(rho_conditions(rdms, 2),:,:))));
rho2z = @(x) sqrt((numel(ind)-3)/1.06) * x; % transform Fisher zr into z for obtaining p, see https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient; Fieller et al., 1957
zmap = rho2z(diffmap);

    for j = 1:length(pn)
    threshed{j} = zmap .*pn(j) >= 1.96;
    [cnum, N] = bwlabel(threshed{j}, 4); % 4-connected (no vertical connection)
    clusters{j}    = arrayfun(@(x) find(cnum==x)', 1:N, 'UniformOutput', false); % sorted labels
    clusterSize{j} = arrayfun(@(x) numel([x{:}]), clusters{j});
    pSigClusters{j} = arrayfun(@(x) sum([mxc{:}] > x) / numel([mxc{:}]), clusterSize{j});
    sigIndex{j} = find(pSigClusters{j} < .025);
    cl = clusters{j};
    si = sigIndex{j};
    if isempty(si)
    clcell{j} = [];
    else
    clcell{j} = [cl{si}]; % continue here
    end
    disp(num2str(pn(j)));
    disp(['Size: ', num2str(clusterSize{j})]);
    disp(['p: ', num2str(pSigClusters{j})]);
    disp('');
    end
    clcellrdms{rdms} = clcell;
    threshedrdms{rdms} = threshed;
    sigTFR = zeros(size(zmap)); sigTFR([clcell{:}]) = zmap([clcell{:}]);

subplot(3, 2, 2*(rdms-1)+1);
imagesc([-0.6:0.02:1.18],30:-1:4, diffmap, [-0.1 0.1]); % vormals [-0.06 0.06]
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Spearman''s rho';
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title(rho_labels{rdms});

subplot(3, 2, 2*(rdms-1)+2);
imagesc([-0.6:0.02:1.18],30:-1:4, sigTFR, [-5 5]); % vormals [-3 3]
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'z score';
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title('Cluster correction');
end

width=1200;
height=800;
set(gcf,'position',[rectFig(1),rectFig(2)-height/2,width,height], 'color', 'white');


intersct = NaN(size(zmap));
[i1, i2]  = deal(zeros(size(zmap)));
i1(clcellrdms{2}{1}) = 1;
i2(clcellrdms{3}{1}) = 2;
figure(5);
imagesc([-0.6:0.02:1.18],30:-1:4, i1+i2, [0 3]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
yline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
title('Intersection of effects');
text(0, 5, 'light blue: basic vs. super', 'Color', 'white');
text(0, 7, 'green: sub vs. super', 'Color', 'white');
text(0, 9, 'yellow: intersection', 'Color', 'white');
width=600;
height=300;
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

%%%
% clst = clusters{1};
% [rw,cl] = ind2sub([size(da_grandmean,3), size(da_grandmean,4)], clst(2));
% pltindex = sortrows([freqaxis(rw); timeaxis(cl)]'); % indices for plotting

% strplt  = sortrows(pltindex, [2 1]);
% strplt(diff(strplt(:,1))==1 & diff(strplt(:,2)) == 0,:) = [];
% x = strplt(:,2); y = abs(34-strplt(:,1)); % correct y because of scaling
% 
% pltxy = [x(1), y(1)-0.5;
%          x(2), y(1)-0.5;
%          x(2), y(2)-0.5;
%          x(7), y(7)-0.5;
%          x(7), y(8)-0.5;
%          x(8), y(8)-0.5;
%          x(8), y(8)+0.5;
%          x(1), y(8)+0.5;
%          x(1), y(8)-0.5;];
% hold on 
% plot(pltxy(:,1), pltxy(:,2), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');


