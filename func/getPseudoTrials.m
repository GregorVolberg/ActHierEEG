function [pTrials] = getPseudoTrials(data, nbins, time, bsltime)
% randomize trials into approx. equally sized bins for pseudotrials;


bindices = getbins(data, nbins);

downsampledTime = time(1:2:181-1);
bsltimebins = [find(downsampledTime == bsltime(1)), find(downsampledTime == bsltime(2))];

pTrials = NaN(nbins, size(data,2), size(data,3), length(downsampledTime));

    for pseudoTrial = 1:numel(bindices) % create pseudotrials
        desc = squeeze(nanmean(data(bindices{pseudoTrial},:,:,:),1)); % 62 x 27 x 181
        % downsample by averaging
        d1 = desc(:,:,1:2:181-1);
        d2 = desc(:,:,2:2:181);
        downsampledDesc = (d1+d2)/2;
        meanbsl = nanmean(downsampledDesc(:,:,bsltimebins(1):bsltimebins(2)),3); % 62 x 27
        descbsl = 10*log10(downsampledDesc ./ meanbsl);
        %descbsl = ft_freqbaseline(cfgbsl, desc);
        pTrials(pseudoTrial,:,:,:) = descbsl;
    end

end
