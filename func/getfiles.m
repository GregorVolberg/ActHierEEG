function [ffiles] = getfiles(logpath, rawpath, cleanpath)
% raw files
info = dir([logpath, '*_01.mat']);
vhdr = dir([rawpath, '*.vhdr']);
ffiles.info = sort(strcat({info.folder}, filesep, {info.name}))'; clear info
ffiles.vhdr = sort(strcat({vhdr.folder}, filesep, {vhdr.name}))'; clear vhdr

% cleaned files
tmp     = dir([cleanpath, '*.mat']);
tmpchar = char({tmp.name});
vps     = cellstr(tmpchar(:,1:5));
cleanfiles = strcat({tmp.folder}, filesep, {tmp.name})'; clear tmp

% reduce raw file list to cleaned file list
matchedFiles = contains(ffiles.vhdr, vps);
ffiles.info(matchedFiles==0) = [];
ffiles.vhdr(matchedFiles==0) = [];
ffiles.clean                 = cleanfiles;
ffiles.vps                   = vps;

% check
%tmp = char(ffiles.info);  ffi = cellstr(tmp(:,54:58));
%tmp = char(ffiles.vhdr);  ffv = cellstr(tmp(:,54:58));
%tmp = char(ffiles.clean); ffc = cellstr(tmp(:,67:71));
%vp = char(ffiles.vps);
%char(ffi)==char(ffv)
%char(ffi)== vp
%etc
end