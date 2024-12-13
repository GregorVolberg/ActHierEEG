%addpath('X:/Volberg/m-lib/fieldtrip')
%addpath('/nas_rz_share/data/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');
addpath(genpath('./func/libsvm'));

% partly adopted from https://github.com/anonymturtle/VCR_infant/tree/main/code
% "Visual category representation in the infant brain"

tfrpath    = './tfr/';
dapath     = './DAtfr2/';

tmp = dir([tfrpath, 'SUB*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

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

% for recoding trial info, according to Tonghe's script
orderOfActions=[49,50,51,52,53,54,13,14,15,16,17,18,37,38,39,40,41,42,55,56,...
    57,58,59,60,7,8,9,10,11,12,61,62,63,64,65,66,1,2,3,4,5,6,43,44,45,46,47,48,...
    19,20,21,22,23,24,25,26,27,28,29,30,67,68,69,70,71,72,31,32,33,34,35,36];

%% subject loop

for vp = 18:numel(files.name)
%vp=1

eeg = load (files.name{vp});
fprintf([files.vps{vp}, ' loaded\n']);
%bsltimebins = [find(eeg.tfr.time == baseline(1)), find(eeg.tfr.time == baseline(2))];
% results matrix pre-allocation
DAmean = nan(nconditions, nconditions, size(eeg.tfr.powspctrm, 3), (size(eeg.tfr.powspctrm, 4)-1)/2);

% recode trial code so that it matches row and columns in RDM 
% corresponding stimulus names see bottom of script
trialCodes = orderOfActions(eeg.tfr.trialinfo(:,4));

for condRun = 1:size(conds,1)
%    fprintf(['\nCondition row ', num2str(condRun)]);       
    trialsA = find(trialCodes == conds(condRun, 1)); % condition A
    selA = eeg.tfr.powspctrm(trialsA, :,:,:);
    trialsB = find(trialCodes == conds(condRun, 2)); % condition B
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
end 

  % 72×1 cell array
  % 
  %   {'./stimuli/Motorrollerfahren1.jpg'}
  %   {'./stimuli/Motorrollerfahren2.jpg'}
  %   {'./stimuli/Motorrollerfahren3.jpg'}
  %   {'./stimuli/Motorrollerfahren4.jpg'}
  %   {'./stimuli/Motorrollerfahren5.jpg'}
  %   {'./stimuli/Motorrollerfahren6.jpg'}
  %   {'./stimuli/Fahrradfahren1.jpg'    }
  %   {'./stimuli/Fahrradfahren2.jpg'    }
  %   {'./stimuli/Fahrradfahren3.jpg'    }
  %   {'./stimuli/Fahrradfahren4.jpg'    }
  %   {'./stimuli/Fahrradfahren5.jpg'    }
  %   {'./stimuli/Fahrradfahren6.jpg'    }
  %   {'./stimuli/Kraulenschwimmen1.jpg' }
  %   {'./stimuli/Kraulenschwimmen2.jpg' }
  %   {'./stimuli/Kraulenschwimmen3.jpg' }
  %   {'./stimuli/Kraulenschwimmen4.jpg' }
  %   {'./stimuli/Kraulenschwimmen5.jpg' }
  %   {'./stimuli/Kraulenschwimmen6.jpg' }
  %   {'./stimuli/Rückenschwimmen1.jpg'  }
  %   {'./stimuli/Rückenschwimmen2.jpg'  }
  %   {'./stimuli/Rückenschwimmen3.jpg'  }
  %   {'./stimuli/Rückenschwimmen4.jpg'  }
  %   {'./stimuli/Rückenschwimmen5.jpg'  }
  %   {'./stimuli/Rückenschwimmen6.jpg'  }
  %   {'./stimuli/Biertrinken1.jpg'      }
  %   {'./stimuli/Biertrinken2.jpg'      }
  %   {'./stimuli/Biertrinken3.jpg'      }
  %   {'./stimuli/Biertrinken4.jpg'      }
  %   {'./stimuli/Biertrinken5.jpg'      }
  %   {'./stimuli/Biertrinken6.jpg'      }
  %   {'./stimuli/Wassertrinken1.jpg'    }
  %   {'./stimuli/Wassertrinken2.jpg'    }
  %   {'./stimuli/Wassertrinken3.jpg'    }
  %   {'./stimuli/Wassertrinken4.jpg'    }
  %   {'./stimuli/Wassertrinken5.jpg'    }
  %   {'./stimuli/Wassertrinken6.jpg'    }
  %   {'./stimuli/Apfelessen1.jpg'       }
  %   {'./stimuli/Apfelessen2.jpg'       }
  %   {'./stimuli/Apfelessen3.jpg'       }
  %   {'./stimuli/Apfelessen4.jpg'       }
  %   {'./stimuli/Apfelessen5.jpg'       }
  %   {'./stimuli/Apfelessen6.jpg'       }
  %   {'./stimuli/Kuchenessen1.jpg'      }
  %   {'./stimuli/Kuchenessen2.jpg'      }
  %   {'./stimuli/Kuchenessen3.jpg'      }
  %   {'./stimuli/Kuchenessen4.jpg'      }
  %   {'./stimuli/Kuchenessen5.jpg'      }
  %   {'./stimuli/Kuchenessen6.jpg'      }
  %   {'./stimuli/Fensterputzen1.jpg'    }
  %   {'./stimuli/Fensterputzen2.jpg'    }
  %   {'./stimuli/Fensterputzen3.jpg'    }
  %   {'./stimuli/Fensterputzen4.jpg'    }
  %   {'./stimuli/Fensterputzen5.jpg'    }
  %   {'./stimuli/Fensterputzen6.jpg'    }
  %   {'./stimuli/Geschirrabwaschen1.jpg'}
  %   {'./stimuli/Geschirrabwaschen2.jpg'}
  %   {'./stimuli/Geschirrabwaschen3.jpg'}
  %   {'./stimuli/Geschirrabwaschen4.jpg'}
  %   {'./stimuli/Geschirrabwaschen5.jpg'}
  %   {'./stimuli/Geschirrabwaschen6.jpg'}
  %   {'./stimuli/Zähneputzen1.jpg'      }
  %   {'./stimuli/Zähneputzen2.jpg'      }
  %   {'./stimuli/Zähneputzen3.jpg'      }
  %   {'./stimuli/Zähneputzen4.jpg'      }
  %   {'./stimuli/Zähneputzen5.jpg'      }
  %   {'./stimuli/Zähneputzen6.jpg'      }
  %   {'./stimuli/Gesichtwaschen1.jpg'   }
  %   {'./stimuli/Gesichtwaschen2.jpg'   }
  %   {'./stimuli/Gesichtwaschen3.jpg'   }
  %   {'./stimuli/Gesichtwaschen4.jpg'   }
  %   {'./stimuli/Gesichtwaschen5.jpg'   }
  %   {'./stimuli/Gesichtwaschen6.jpg'   }
