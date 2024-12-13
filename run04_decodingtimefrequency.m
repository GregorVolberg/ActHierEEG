%addpath('X:/Volberg/m-lib/fieldtrip')
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;
tic
addpath('./func');
addpath(genpath('./func/libsvm'));

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

tfrpath    = './tfr/';
dapath     = './DA/';

tmp = dir([tfrpath, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

%% Download dataset (if necessary) and add VCR_infant to the MATLAB path
%setup([3,4]);


%% Sort trials into conditions
% 72 pictures; each one condition; sorted so that 1:6 belong to action1,
% 7:12 action2 etc

% how many pseudotrials?
nbins = 4;

baseline           = [-0.5 -0.1];

%cfgbsl.baselinetype       = 'db';

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

for vp = 18:numel(files.name)
%vp=1

eeg = load (files.name{vp});
fprintf([files.vps{vp}, ' loaded\n']);
%bsltimebins = [find(eeg.tfr.time == baseline(1)), find(eeg.tfr.time == baseline(2))];
% results matrix pre-allocation
DAmean = nan(nconditions, nconditions, size(eeg.tfr.powspctrm, 3), (size(eeg.tfr.powspctrm, 4)-1)/2);

for condRun = 1:size(conds,1)
%    fprintf(['\nCondition row ', num2str(condRun)]);       
    trialsA = find(eeg.tfr.trialinfo(:,4) == conds(condRun, 1)); % condition A
    selA = eeg.tfr.powspctrm(trialsA, :,:,:);
    trialsB = find(eeg.tfr.trialinfo(:,4) == conds(condRun, 2)); % condition B
    selB = eeg.tfr.powspctrm(trialsB, :,:,:); 

    DA  = zeros(size(eeg.tfr.powspctrm, 3), (size(eeg.tfr.powspctrm, 4)-1)/2); % pre-allocate; decimate time points
    %clear eeg
    
    for permrun = 1:numruns

    pseudoA = getPseudoTrials(selA, nbins, eeg.tfr.time, baseline);
    pseudoB = getPseudoTrials(selB, nbins, eeg.tfr.time, baseline);
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

toc

end 
% 
% % Check dataCell
% %dataCell
% %%
% 
% %% Time-frequency decomposition
% %[dBPowerCell,convo] = timefrexdecomp(dataCell);
% %clear dataCell;
% 
% % Check dBPowerCell and convo
% %dBPowerCell
% %convo
% %%
% 
% %% Define time points and frequencies of interest
% 
% % Time points of interest
% %TOI = dsearchn(convo.times',[-200, 1000]');
% %dBPowerCell = cellfun(@(x) x(:,:,TOI(1):TOI(2)), dBPowerCell,'UniformOutput', 0);
% 
% % Frequencies of interest
% %FOI = convo.mor_frex;
% %frequencyF = length(FOI);
% %%
% 
% %% Decode object categories in time and frequency
% 
% % Time- and frequency- resolved classification analysis
% % The procedure is equivalent to the time-resolved classification analysis
% % (see example_decodingintime.m) with the only difference that
% % classification is conducted on power value patterns instead of 
% % raw activation value patterns. 
% % The analysis is conducted for each frequency bin separately. 
% 
% % Pre-allocate cell array
% %DA_frexCell = cell(frequencyF,1);
% 
% for frexF = 1:frequencyF %Loop through frequencies
%     
%     % Define dataCell for decoding, i.e., power data in one frequency step
%  %   dataCell = {dBPowerCell{1:4,frexF}};
%  %   dataCell = reshape(dataCell,[4,1]);
%     
%     % For time and memory efficiency, we downsampled dB power data to a 
%     % sampling rate of 50 Hz by averaging the raw EEG data in 20ms bins.
%  %   origTempRes = 500; % 500 Hz
%  %   resTempRes = 50; % 50 Hz
%  %   timePoints2Average = origTempRes/resTempRes;
%  %   [dataCell, timepointT] = averagetimepoints(dataCell, timePoints2Average);
%     % timepointT: Time points of EEG epoch (from -200ms to +1000ms) in 20ms-steps
%     
%     % Decoding is done pair-wise for all pairs of conditions.
%     % For M conditions, you will have to do ((M*(M-1))/2 (i.e., possible
%     % pair-wise condition combinations) condition-pair specific classifications.
%     % The result you could get for a pairwise classification would be 0, 50 or
%     % 100% correct with 50% being chance level.
%     
%     % NOTE: Randomization is used within averagetrials.m so results might
%     % differ between computations.
% 
%     % We repeat the training & testing procedure "permutationX" times average
%     % the accuracies of all repetitions to get the mean decoding accuracy.
%     permutationX = 10;
%     % Theoretically, the more repetitions, the finer the sampling of a real
%     % state of things. 100 permutations is usually enough, one can however
%     % determine empirically the number for a given dataset.
%     % Here, we use 10 to speed up the example.
%     
%     % Object categories as conditions
%     conditionM = 4;
%     
%     % Pre-allocate result matrix
%     DA = nan(permutationX, conditionM, conditionM, timepointT);
%     % This analysis yields a M*M decoding accuracy matrix indexed
%     % in rows and columns by the conditions compared for all time point
%     % combinations from -200ms to +1,000ms.
%     
%     for permX = 1:permutationX % Loop through repetitions
%         
%         % To increase the signal-to-noise ratio (SNR),
%         % we randomly assigned raw trials into N bins of approximately
%         % equal size each and averaged them into pseudo-trials. 
%         % In this example, N = 4.
%         pseudoTrialN = 4;
%         pseudoData = averagetrials(dataCell, pseudoTrialN);
%         
%         for condA = 1:conditionM % Loop for condition A
%             for condB = condA+1:conditionM % Loop for condition B
%                 for timeT = 1:timepointT % Loop for time point A
%                     
%                     % Implement leave-one-pseudo-trial-out cross validated
%                     % classification approach.
%                     
%                     % We trained the SVM classifier to perform pair-wise decoding on any 
%                     % two conditions using three of the four pseudo-trials 
%                     % for training at time point A.
%                     training_data = ...
%                         double([squeeze(pseudoData(condA, 1:end - 1, :, timeT)); ...
%                         squeeze(pseudoData(condB, 1:end - 1, :, timeT))]);
%                     labels_train = [ones(pseudoTrialN-1,1);...
%                         2*ones(pseudoTrialN - 1,1)];
%                     
%                     % Train model
%                     model = svmtrain(labels_train, training_data, '-s 0 -t 0 -q');
%                     
%                     % We used the fourth left-out pseudo-trial
%                     % at time point B for testing, yielding classification
%                     % accuracy (chance level 50%) as a result.
%                     testing_data = double([squeeze(pseudoData(condA, end, :, timeT))'; ...
%                         squeeze(pseudoData(condB, end, :, timeT))']);
%                     labels_test = [1;2];
%                     
%                     % Test model
%                     [~, accuracy, ~] = svmpredict(labels_test, testing_data, model, '-q');
%                     
%                     % Assign the decoding result into DA matrix
%                     DA(permX, condB, condA, timeT) = accuracy(1);
%                     
%                 end
%             end
%         end
%     end
%     
%     % Average the DA matrix across repetitions
%     DA_frexCell{frexF,1} = nanmean(DA,1);
% end
% 
% DA_timefrex = cat(1,DA_frexCell{:}); 
% % 4D matrix: frequency x condition x condition x time point
% 
% % Display run time
% disp("Decoding done.")
% runTime_minutes = toc/60
% %%
% 
% %% Plot EEG dB power after decomposition
% 
% % We average dB power values across channels, trials
% meandBPowerCell = cellfun(@(x) mean(mean(x,1),2),dBPowerCell,'UniformOutput',0);
% meandBPowerMat = cell2mat(meandBPowerCell);
% 
% figure(1);
% imagesc(-200:2:1000,1:30,squeeze(nanmean(meandBPowerMat,1)));
% set(gca,'YTick',[1,6,11,16,21,26,30], ...
%     'YTickLabel', round(FOI([1,6,11,16,21,26,30]),2));
% axis xy; axis tight;
% xlim([-200,1000]);
% xlabel('Time (ms)'); ylabel('Frequency (Hz)');
% CH = colorbar('eastoutside'); CH.Label.String = 'Power (dB)';
% title('Decomposition: average EEG signal power');
% % Adjust the position and size of figure
% rectFig = get(gcf,'position');
% width=600;
% height=300;
% set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
% %%
% 
% 
% %% Plot decoding matrix at (100ms, 2Hz)
% 
% % For each time-frequency combination (e.g., tx = 100ms; fy = 2Hz),
% % we get a decoding accuracy matrix of size 4 × 4, with rows and columns
% % indexed by the conditions classified.
% % The matrix is symmetric across the diagonal, with the diagonal itself being undefined.
% % This procedure yielded one decoding matrix for every time point.
% tx = dsearchn((-200:20:1000)', 100);
% fy = dsearchn(FOI', 2);
% DA_mean_matrix = squeeze(DA_timefrex(fy,:,:,tx));
% DecodingMatrix = triu(DA_mean_matrix.',1) + tril(DA_mean_matrix);
% 
% figure(2);
% imagesc(1:4,1:4,DecodingMatrix);
% set(gca, 'xtick',1:4, 'xticklabel',{'toy','body','house','face'},...
%     'ytick',1:4, 'yticklabel',{'toy','body','house','face'});
% axis equal; axis tight;
% xlabel('Condition'); ylabel ('Condition');
% CH = colorbar('southoutside');
% CH.Label.String = 'Decoding accuracy (%)';
% title('Decoding matrix at (100ms,2Hz)');
% % Adjust the position and size of figure
% rectFig = get(gcf,'position');
% width=600;
% height=300;
% set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
% %%
% 
% %% Plot time-frequency decoding matrix
% figure(3); 
% imagesc(-200:20:1000, 1:30, squeeze(nanmean(nanmean(DA_timefrex,2),3)));
% set(gca,'YTick',[1,6,11,16,21,26,30], ...
%     'YTickLabel', round(FOI([1,6,11,16,21,26,30]),2));
% axis xy; axis tight;
% xlim([-200,1000]);
% xlabel('Time (ms)'); ylabel('Frequency (Hz)');
% CH = colorbar('eastoutside');
% CH.Label.String = 'Decoding accuracy (%)';
% title('Time-frequency decoding matrix');
% % Adjust the position and size of figure
% rectFig = get(gcf,'position');
% width=600;
% height=300;
% set(gcf,'position',[rectFig(1),rectFig(2),width,height], 'color', 'white');
%%
