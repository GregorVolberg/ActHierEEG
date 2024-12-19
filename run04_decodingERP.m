% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

addpath('./func');
addpath(genpath('./func/libsvm'));

erppath    = '../data/erp/';
dapath     = '../data/DAerp/';
files      = get_erpfiles(erppath);

%% Sort trials into conditions
% 72 pictures; each one condition; sorted so that 1:6 belong to action1,
% 7:12 action2 etc

% matrix of indices for which DA is to be computed
nconditions = 72;
conds = get_rc_indices(nconditions);

% randomization parameters
nbins = 4; % how many pseudotrials?
numruns = 100; % runs for testing
rng(12); % seed for random number

%% participant loop
for vp = 21:numel(files.name)

erp = importdata(files.name{vp});
fprintf([files.vps{vp}, ' loaded, trialinfo OK\n']);

% results matrix pre-allocation
DA = nan(nconditions, nconditions, size(erp.trial, 3)); % 72 x 72 x 601

for condRun = 1:size(conds,1)
    trialsA = find(erp.trialinfo(:,5) == conds(condRun, 1)); % condition A
    selA = erp.trial(trialsA, :,:);
    trialsB = find(erp.trialinfo(:,5) == conds(condRun, 2)); % condition B
    selB = erp.trial(trialsB, :,:); 

    DArun  = zeros(1, size(erp.trial, 3)); % pre-allocate 1 x 601
    
    for permrun = 1:numruns

    pseudoA = getPseudoTrialsERP(selA, nbins);
    pseudoB = getPseudoTrialsERP(selB, nbins);
    tmpDA   = NaN(size(DArun)); %pre-allocate
    
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
    DArun = DArun + tmpDA;
    end %permrun
    
    DA(conds(condRun, 1), conds(condRun, 2), :) = DArun./numruns;
end %condRun


save([dapath, files.vps{vp}, '_DA.mat'], 'DA');

end 
