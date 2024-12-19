% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

% paths and folders
addpath('./func');
addpath(genpath('./func/libsvm'));

tfrpath    = '../data/tfr/';
dapath     = '../data/DAtfr/';
files      = get_tfrfiles(tfrpath);

% matrix of indices for which DA is to be computed
nconditions = 72;
conds = get_rc_indices(nconditions);

% randomization parameters
nbins = 4; % how many pseudotrials?
numruns = 100; % for testing
rng(12); % seed for random number

% baseline specification; used in getPseudoTrialsTFR()
baseline           = [-0.5 -0.1]; 

%% subject loop
for vp = 1:numel(files.name)
eeg = load (files.name{vp});
fprintf([files.vps{vp}, ' loaded\n']);

DA = nan(nconditions, nconditions, size(eeg.tfr.powspctrm, 3), (size(eeg.tfr.powspctrm, 4)-1)/2);

for condRun = 1:size(conds,1)
    trialsA = find(eeg.tfr.trialinfo(:,4) == conds(condRun, 1)); % condition A
    selA = eeg.tfr.powspctrm(trialsA, :,:,:);
    trialsB = find(eeg.tfr.trialinfo(:,4) == conds(condRun, 2)); % condition B
    selB = eeg.tfr.powspctrm(trialsB, :,:,:); 

    DArun  = zeros(size(eeg.tfr.powspctrm, 3), (size(eeg.tfr.powspctrm, 4)-1)/2); % pre-allocate; decimate time points
        
    for permrun = 1:numruns
        pseudoA = getPseudoTrialsTFR(selA, nbins, eeg.tfr.time, baseline); % 
        pseudoB = getPseudoTrialsTFR(selB, nbins, eeg.tfr.time, baseline);
        tmpDA  = NaN(size(DArun)); %pre-allocate
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
                tmpDA(freqbin, timeP) = accuracy(1);
            end %time
        end %freq
    DArun = DArun + tmpDA;
    end %permrun
    DA(conds(condRun, 1), conds(condRun, 2), :, :) = DArun./numruns;
end %condRun

save([dapath, files.vps{vp}, '_DA.mat'], 'DA');
end 
