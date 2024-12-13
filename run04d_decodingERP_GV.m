% same data as fpr TFR analysis

addpath('./func');
addpath(genpath('./func/libsvm'));

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

erppath    = './erp/';
dapath     = './DAerp2/';

tmp = dir([erppath, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

%% Sort trials into conditions
% 72 pictures; each one condition; sorted so that 1:6 belong to action1,
% 7:12 action2 etc

% how many pseudotrials?
nbins = 4;

baseline           = [-0.2 0]; % already subtracted in single trials

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

for vp = 1:numel(files.name)
eeg = load (files.name{vp});

% results matrix pre-allocation
DAmean = nan(nconditions, nconditions, size(eeg.erp.trial, 3));

for condRun = 1:size(conds,1)
    trialsA = find(eeg.erp.trialinfo(:,4) == conds(condRun, 1)); % condition A
    selA = eeg.erp.trial(trialsA, :,:);
    trialsB = find(eeg.erp.trialinfo(:,4) == conds(condRun, 2)); % condition B
    selB = eeg.erp.trial(trialsB, :,:); 

    DA  = zeros(1, size(eeg.erp.trial, 3)); % pre-allocate; decimate time points
    %clear eeg
    
    for permrun = 1:numruns

    pseudoA = getPseudoTrialsERP(selA, nbins);
    pseudoB = getPseudoTrialsERP(selB, nbins);
    tmpDA  = NaN(size(DA)); %pre-allocate
    
         for timeP = 1:size(pseudoA, 3)
            
        trainingdata = double(squeeze([pseudoA(1:nbins-1, :, timeP);
                        squeeze(pseudoB(1:nbins-1, :, timeP))]));
        labels_training = [ones(nbins-1,1);...
                            2*ones(nbins-1,1)];                
        model = svmtrain(labels_training, trainingdata, '-s 0 -t 0 -q');

        testingdata = double([squeeze(pseudoA(end, :, timeP)); ...
                            squeeze(pseudoB(end, :, timeP))]);
        labels_testing = [1; 2];
        [~, accuracy, ~] = svmpredict(labels_testing, testingdata, model, '-q');

         tmpDA(timeP) = accuracy(1);
        end %time
 
    DA = DA + tmpDA;
    end %permrun
    
    DAmean(conds(condRun, 1), conds(condRun, 2), :, :) = DA./numruns;
end %condRun


save([dapath, files.vps{vp}, 'DAmean.mat'], 'DAmean');

end 

