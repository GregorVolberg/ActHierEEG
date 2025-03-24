gvpath = './DAerp/';
tzpath = '../Results/RSA/EEG_based_RDM/';

tmp = dir([gvpath, '*DAmean.mat']);
tmpchar = char({tmp.name});
gvfiles.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
gvfiles.vps  = cellstr(tmpchar(:,1:5));

tmp = dir([tzpath, '*_decAll.mat']);
tmpchar = char({tmp.name});
tzfiles.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
tzfiles.vps  = cellstr(tmpchar(:,1:5));

% same participants in same order
all(all(char(tzfiles.vps) == char(gvfiles.vps)))

[gv, tz] = deal(nan(22, 72, 72, 600));
for vp = 1:numel(gvfiles.vps)
    gv(vp, :,:,:) = importdata(gvfiles.name{vp});
    tmp = load(tzfiles.name{vp});
    tz(vp, :,:,:) = tmp.res.diss; 
end

nconditions = size(gv,2);
ind = find(tril(ones(nconditions, nconditions), -1) == 1);

[DAgv, DAtz] = deal(nan(22, 600));
for subject = 1 : size(gv, 1)
    for timepoint = 1: size(gv, 4)
        tmp  = squeeze(gv(subject, :, :, timepoint));
        DAgv(subject, timepoint) = mean(tmp(ind));
        tmp  = squeeze(tz(subject, :, :, timepoint));
        DAtz(subject, timepoint) = mean(tmp(ind));
    end
end

figure(1);
plot(DAgv'); hold on; plot(squeeze(mean(DAgv, 1))', 'Linewidth', 2, 'Color', 'black');
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([40 80]); xlabel('Time (s)'); ylabel('Decoding accuracy Gregor');
figure(2);
plot(DAtz'); hold on; plot(squeeze(mean(DAtz, 1))', 'Linewidth', 2, 'Color', 'black');
set(gca,'XTick',[1:100:600], ...
     'XTickLabel', [-0.2:0.2:1]);%{'subordinate', 'basic', 'superordinate'}, ...
ylim([0.4 0.8]); xlabel('Time (s)'); ylabel('Decoding accuracy Tonghe');

