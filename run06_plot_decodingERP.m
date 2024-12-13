%addpath('./func');

%erppath     = './erp/';
dapath      = './DAerp/';
%dbpowfile   = [erppath, 'meandbpow.mat'];
dafile      = [dapath, 'DA_all.mat'];
permfile    = [dapath, 'permttest.mat'];
clusterfile = [dapath, 'clustermat.mat'];

% Plot DA
load(dafile);
meanDA  = squeeze(mean(DA_all, 1));
figure(1);
plot([-0.2:0.002:0.998], meanDA');
ylim([45, 65]);
set(gca,'YTick',[50:5:65], ...
    'YTickLabel', [50:5:65]);
xlabel('Time (s)'); ylabel('Decoding Accuracy (%)');
title('Average decoding accuracy, ERP a la Tonghe');
rectFig = get(gcf,'position');
width=600;
height=300;
yline(50, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
 
% plot stats
%load(permfile);
load(clusterfile);
hold on
line(clustermat.time, [49 49], color = 'red')

