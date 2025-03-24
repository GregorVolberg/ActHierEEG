function [d_accuracy] = get_decodingaccuracy_fromfiles(da_path)

tmp = dir([da_path, '*DAmean.mat']);
tmpchar = char({tmp.name});
dafiles.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
dafiles.vps  = cellstr(tmpchar(:,1:5));

for vp = 1:numel(dafiles.vps)
    load(dafiles.name{vp});
    da(vp, :,:,:,:) = DAmean;
end

d_accuracy.files    = dafiles.name;
d_accuracy.subjects = dafiles.vps;
d_accuracy.data     = da;
end
