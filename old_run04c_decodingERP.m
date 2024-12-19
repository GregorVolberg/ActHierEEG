%addpath('X:/Volberg/m-lib/fieldtrip')
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;

addpath('./func');
addpath(genpath('./func/libsvm'));

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

erppath1    = '../Results/Preprocessing/';
erppath2    = './cleaned/';
dapath     = './DAerp3/';

tmp = dir([erppath1, '*.mat']);
tmpchar = char({tmp.name});
files.name1 = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
tmp = dir([erppath2, '*.mat']);
tmpchar = char({tmp.name});
files.name2 = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
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

% for recoding trial info, according to Tonghe's script
orderOfActions=[49,50,51,52,53,54,13,14,15,16,17,18,37,38,39,40,41,42,55,56,...
    57,58,59,60,7,8,9,10,11,12,61,62,63,64,65,66,1,2,3,4,5,6,43,44,45,46,47,48,...
    19,20,21,22,23,24,25,26,27,28,29,30,67,68,69,70,71,72,31,32,33,34,35,36];

%% subject loop

for vp = 1:numel(files.name1)
%vp=1

eeg = load (files.name1{vp});
cln = load (files.name2{vp});
if all(eeg.timelock.trialinfo == cln.data_ica_cleaned.trialinfo(:,1))
    eeg.timelock.trialinfo = cln.data_ica_cleaned.trialinfo;
else
    error("trialinfo mismatch");
end

fprintf([files.vps{vp}, ' loaded, trialinfo OK\n']);

% results matrix pre-allocation
DAmean = nan(nconditions, nconditions, size(eeg.timelock.trial, 3));

% recode trial code so that it matches row and columns in RDM

trialCodes = orderOfActions(eeg.timelock.trialinfo(:,4));

for condRun = 1:size(conds,1)
%    trialsA = find(eeg.timelock.trialinfo(:,4) == conds(condRun, 1)); % condition A
    trialsA = find(trialCodes == conds(condRun, 1)); % condition A
    selA = eeg.timelock.trial(trialsA, :,:);
%    trialsB = find(eeg.timelock.trialinfo(:,4) == conds(condRun, 2)); % condition B
    trialsB = find(trialCodes == conds(condRun, 2)); % condition B
    selB = eeg.timelock.trial(trialsB, :,:); 

% % match trial count per condition
% if size(selA,1) > size(selB,1)
%     tmpsel         = randperm(size(selA,1));
%     selA = selA(tmpsel(1:size(selB,1)), :, :, :);
% elseif size(selB,1) > size(selA,1)
%     tmpsel         = randperm(size(selB,1));
%     selB = selB(tmpsel(1:size(selA,1)), :, :, :);
% end

    DA  = zeros(1, size(eeg.timelock.trial, 3)); % pre-allocate; decimate time points
    %clear eeg
    
    for permrun = 1:numruns

    pseudoA = getPseudoTrialsERP(selA, nbins);
    pseudoB = getPseudoTrialsERP(selB, nbins);
    tmpDA  = NaN(size(DA)); %pre-allocate
    
    %for freqbin = 1:size(pseudoA, 3)
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

        %DA(permrun, condA, condB, freqbin, timeP) = accuracy(1);
        tmpDA(timeP) = accuracy(1);
        end %time
    %end %freq
    DA = DA + tmpDA;
    end %permrun
    
    DAmean(conds(condRun, 1), conds(condRun, 2), :, :) = DA./numruns;
end %condRun


save([dapath, files.vps{vp}, 'DAmean.mat'], 'DAmean');

end 
