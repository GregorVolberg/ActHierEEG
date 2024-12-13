%addpath('X:/Volberg/m-lib/fieldtrip')
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;
tic
addpath('./func');
addpath(genpath('./func/libsvm'));

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

tfrpath    = './itc/';
dapath     = './DAitc/';

tmp = dir([tfrpath, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

%% Sort trials into conditions
% 72 pictures; each one condition; sorted so that 1:6 belong to action1,
% 7:12 action2 etc

% how many pseudotrials?
nbins = 4;

baseline           = [-0.5 -0.1];

% matrix of indices for which DA is to be computed
nconditions = 72;
mat = tril(ones(nconditions, nconditions), -1);
ind = find(mat == 1);
[r, c] = ind2sub(size(mat), ind);
conds = [r, c]; 

% randomization parameters
numruns = 100; % for testing
rng(12); % seed for random number


%% subject loop

for vp = 20:numel(files.name)
%vp=1

eeg = load (files.name{vp});
fprintf([files.vps{vp}, ' loaded\n']);
%bsltimebins = [find(eeg.tfr.time == baseline(1)), find(eeg.tfr.time == baseline(2))];
% results matrix pre-allocation
DAmean = nan(nconditions, nconditions, size(eeg.itc.fourierspctrm, 3), (size(eeg.itc.fourierspctrm, 4)-1)/2);

for condRun = 1:size(conds,1)
%    fprintf(['\nCondition row ', num2str(condRun)]);       
    trialsA = find(eeg.itc.trialinfo(:,4) == conds(condRun, 1)); % condition A
    selA = eeg.itc.fourierspctrm(trialsA, :,:,:);
    trialsB = find(eeg.itc.trialinfo(:,4) == conds(condRun, 2)); % condition B
    selB = eeg.itc.fourierspctrm(trialsB, :,:,:); 

% match trial count per condition
if size(selA,1) > size(selB,1)
    tmpsel         = randperm(size(selA,1));
    selA = selA(tmpsel(1:size(selB,1)), :, :, :);
elseif size(selB,1) > size(selA,1)
    tmpsel         = randperm(size(selB,1));
    selB = selB(tmpsel(1:size(selA,1)), :, :, :);
end

    DA  = zeros(size(eeg.itc.fourierspctrm, 3), (size(eeg.itc.fourierspctrm, 4)-1)/2); % pre-allocate; decimate time points
    %clear eeg
    
    for permrun = 1:numruns

    pseudoA = getPseudoTrialsITC(selA, nbins, eeg.itc.time, baseline);
    pseudoB = getPseudoTrialsITC(selB, nbins, eeg.itc.time, baseline);
    tmpDA  = NaN(size(DA)); %pre-allocate
    
    for freqbin = 1:size(pseudoA, 3)
        for timeP = 1:size(pseudoA, 4)
            
        trainingdata = double(squeeze([pseudoA(1:nbins-1, :, freqbin, timeP);
                        squeeze(pseudoB(1:nbins-1, :, freqbin, timeP))]));
        labels_training = [ones(nbins-1,1);...
                            2*ones(nbins-1,1)];                
        model = svmtrain(labels_training, trainingdata, '-s 0 -t 0 -q');

        testingdata = double([squeeze(pseudoA(end, :, freqbin, timeP)); ...
                            squeeze(pseudoB(end, :, freqbin, timeP))]);
        labels_testing = [1; 2];
        [~, accuracy, ~] = svmpredict(labels_testing, testingdata, model, '-q');

        %DA(permrun, condA, condB, freqbin, timeP) = accuracy(1);
        tmpDA(freqbin, timeP) = accuracy(1);
        end %time
    end %freq
    DA = DA + tmpDA;
    end %permrun
    
    DAmean(conds(condRun, 1), conds(condRun, 2), :, :) = DA./numruns;
end %condRun


save([dapath, files.vps{vp}, 'DAmean.mat'], 'DAmean');

end 
