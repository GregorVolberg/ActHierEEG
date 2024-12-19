function [files] = get_erpfiles(erppath)
tmp = dir([erppath, 'SUB*erp.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
end