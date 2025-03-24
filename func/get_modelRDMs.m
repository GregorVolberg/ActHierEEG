function [RDM] = get_modelRDMs()
RDMtemplate = ones(72, 72);
RDMsub = RDMtemplate;
k = 0:6:72;
for rn = 1:numel(k)-1
RDMsub(k(rn)+1:k(rn+1), k(rn)+1:k(rn+1)) = 0;
end
%imagesc(RDMsub)

RDMbasic = RDMtemplate;
k = 0:12:72;
for rn = 1:numel(k)-1
RDMbasic(k(rn)+1:k(rn+1), k(rn)+1:k(rn+1)) = 0;
end
%imagesc(RDMbasic)

RDMsuper = RDMtemplate;
k = 0:24:72;
for rn = 1:numel(k)-1
RDMsuper(k(rn)+1:k(rn+1), k(rn)+1:k(rn+1)) = 0;
end
%imagesc(RDMsuper)

RDM = {RDMsub, RDMbasic, RDMsuper};
end
