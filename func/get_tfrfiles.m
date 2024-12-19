function [files] = get_tfrfiles(tfrpath)
tmp = dir([tfrpath, 'SUB*tfr.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
end