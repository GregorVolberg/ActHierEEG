% addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;
%addpath('./func');

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"
tfrpath     = './tfr/';
dapath      = './DA/';
dbpowfile   = [tfrpath, 'meandbpow.mat'];
dafile      = [dapath, 'DA_all.mat'];
permfile    = [dapath, 'permttest.mat'];
clusterfile = [dapath, 'clustermat.mat'];

% Plot EEG dB power
load(dbpowfile);
figure(1);
imagesc(-0.6:0.1:1.2,4:30,flipud(meandbpow), [-1.6,1.6]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Power (dB)';
title('Average EEG power');
rectFig = get(gcf,'position');
width=600;
height=300;
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');

% Plot DA
load(dafile);
meanDA  = squeeze(mean(DA_all, 1));
figure(2);
imagesc([-0.6:0.02:1.18],4:30,flipud(meanDA), [50,58]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside'); CH.Label.String = 'Accuracy (%)';
title('Average decoding accuracy');
rectFig = get(gcf,'position');
width=600;
height=300;
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
 
% plot stats
load(permfile);
load(clusterfile);
figure(3);
imagesc([-0.6:0.02:1.18],30:-1:4,1-permttest.pval, [0.975,1]);
set(gca,'YTick',[4:5:30], ...
    'YTickLabel', [30:-5:5]);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
CH = colorbar('eastoutside');
CH.Label.String = 'p value';
CH.Ticks = [0.975:0.005:1]; 
CH.TickLabels= [0.025:-0.005:0];
title('Significant bins for RSA');
rectFig = get(gcf,'position');
width=600;
height=300;
xline(0, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
hold on

tvec = -0.6:0.02:1.18;
fvec = 4:1:30;
c1 = clustermat(1).cluster;
xx = arrayfun(@(x) tvec([min(c1(c1(:,1)==x,2)), max(c1(c1(:,1)==x,2))])+[-0.01 0.01] , 1:2, 'UniformOutput', false);
x = [xx{:}]';
pltxy = [x(1), 30.5;
         x(2), 30.5;
         x(2), 29.5;
         x(4), 29.5;
         x(4), 28.5;
         x(3), 28.5;
         x(3), 29.5;
         x(1), 29.5;
         x(1), 30.5];

plot(pltxy(:,1), pltxy(:,2), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');

c2 = clustermat(2).cluster;
xx = arrayfun(@(x) tvec([min(c2(c2(:,1)==x,2)), max(c2(c2(:,1)==x,2))])+[-0.01 0.01] , 4:6, 'UniformOutput', false);
x = [xx{:}]';
pltxy = [x(1), 30.5-3;
         x(2), 30.5-3;
         x(2), 29.5-3;
         x(4), 29.5-3;
         x(4), 28.5-3;
         x(6), 28.5-3;
         x(6), 27.5-3;
         x(5), 27.5-3;
         x(5), 28.5-3;
         x(3), 28.5-3;
         x(3), 29.5-3;
         x(1), 29.5-3;
         x(1), 30.5-3];

plot(pltxy(:,1), pltxy(:,2), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
