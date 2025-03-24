function [pTrials] = getPseudoTrialsERP(data, nbins)
% randomize trials into approx. equally sized bins for pseudotrials;


bindices = getbins(data, nbins);

%downsampledTime = time(1:2:181-1);
%bsltimebins = dsearchn(time', bsltime');

pTrials = NaN(nbins, size(data,2), size(data,3));

    for pseudoTrial = 1:numel(bindices) % create pseudotrials
        desc = squeeze(nanmean(data(bindices{pseudoTrial},:,:),1)); % 62 x 600
 %       meanbsl = nanmean(desc(:,bsltimebins(1):bsltimebins(2)),2); % 62 x 1
 %       descbsl = desc - meanbsl;
        %descbsl = ft_freqbaseline(cfgbsl, desc);
        pTrials(pseudoTrial,:,:) = desc;
    end

end
