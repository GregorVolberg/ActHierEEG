function [binindx] = getbins(data, nbins)
% randomize trials into approx. equally sized bins for pseudotrials;
    
%    ntrials = size(data.powspctrm,1);
    ntrials = size(data,1);
    rnd = randperm(ntrials);
    inc = ceil(ntrials/nbins);
    [binindx{1:4}] = deal(rnd(1:inc), ...
                              rnd(inc + 1 : inc * 2), ...
                              rnd(inc * 2 + 1 : inc * 3), ...
                              rnd(inc * 3 + 1 : end));
     binindx = binindx(randperm(numel(binindx))); % so that smallest sized bin is not always last bin

end
