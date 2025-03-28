addpath('X:/Volberg/m-lib/fieldtrip'); ft_defaults;
addpath('./func');

% paths and file names
cleanpath  = '../data/cleaned/';
tfrpath    = '../data/tfr/';

tmp = dir([cleanpath, '*.mat']);
tmpchar = char({tmp.name});
files.name = sort(strcat({tmp.folder}, filesep, {tmp.name}))';
files.vps  = cellstr(tmpchar(:,1:5));
clear tmp tmpchar

% tfr
cfgtfr = [];
cfgtfr.output             = 'pow';
cfgtfr.method             = 'mtmconvol';
cfgtfr.taper              = 'hanning';
cfgtfr.foi                = 4:1:30; % 4 to 40 Hz
cfgtfr.t_ftimwin          = 5./cfgtfr.foi;
cfgtfr.toi                = -0.6:0.01:1.2;%
cfgtfr.pad                = 'nextpow2';
cfgtfr.keeptrials         = 'yes';
% perform baseline correction later on pseudo-trials

% subject loop
for vp = 1:numel(files.name) 
eeg = load (files.name{vp});
tfr =  ft_freqanalysis(cfgtfr, eeg.data_ica_cleaned);
tfr.trialinfo = [tfr.trialinfo, EEG_reshape_to_fMRI_vector(tfr.trialinfo(:,4))];
save([tfrpath, files.vps{vp}, 'tfr.mat'], 'tfr');
clear eeg 
end

% condition matrix in *.trialinfo
% col1:= condition code, taken from Tonghe's segmented EEG data
% col2:= condition code, taken from protocol file (as consistency check)
%        Numbers 1 to 12 refer to different basic level actions
% col3:= pre-stimulus ISI in s, computed from protocol file. Was randomized 1 - 1.5 s
% col4:= number 1 to 72 refering to target stimulus (see protocol files, e. g. ../Data/log/SUB01_01.mat, ExpInfo.stimNames) 
% col5:= re-sorted version of column 4.
%        Numbers 1 to 72 match rows and columns of model RDMs. A list of
%        stimuli is given below:

%     {'./stimuli/Motorrollerfahren1.jpg'}
%     {'./stimuli/Motorrollerfahren2.jpg'}
%     {'./stimuli/Motorrollerfahren3.jpg'}
%     {'./stimuli/Motorrollerfahren4.jpg'}
%     {'./stimuli/Motorrollerfahren5.jpg'}
%     {'./stimuli/Motorrollerfahren6.jpg'}
%     {'./stimuli/Fahrradfahren1.jpg'    }
%     {'./stimuli/Fahrradfahren2.jpg'    }
%     {'./stimuli/Fahrradfahren3.jpg'    }
%     {'./stimuli/Fahrradfahren4.jpg'    }
%     {'./stimuli/Fahrradfahren5.jpg'    }
%     {'./stimuli/Fahrradfahren6.jpg'    }
%     {'./stimuli/Kraulenschwimmen1.jpg' }
%     {'./stimuli/Kraulenschwimmen2.jpg' }
%     {'./stimuli/Kraulenschwimmen3.jpg' }
%     {'./stimuli/Kraulenschwimmen4.jpg' }
%     {'./stimuli/Kraulenschwimmen5.jpg' }
%     {'./stimuli/Kraulenschwimmen6.jpg' }
%     {'./stimuli/Rückenschwimmen1.jpg'  }
%     {'./stimuli/Rückenschwimmen2.jpg'  }
%     {'./stimuli/Rückenschwimmen3.jpg'  }
%     {'./stimuli/Rückenschwimmen4.jpg'  }
%     {'./stimuli/Rückenschwimmen5.jpg'  }
%     {'./stimuli/Rückenschwimmen6.jpg'  }
%     {'./stimuli/Biertrinken1.jpg'      }
%     {'./stimuli/Biertrinken2.jpg'      }
%     {'./stimuli/Biertrinken3.jpg'      }
%     {'./stimuli/Biertrinken4.jpg'      }
%     {'./stimuli/Biertrinken5.jpg'      }
%     {'./stimuli/Biertrinken6.jpg'      }
%     {'./stimuli/Wassertrinken1.jpg'    }
%     {'./stimuli/Wassertrinken2.jpg'    }
%     {'./stimuli/Wassertrinken3.jpg'    }
%     {'./stimuli/Wassertrinken4.jpg'    }
%     {'./stimuli/Wassertrinken5.jpg'    }
%     {'./stimuli/Wassertrinken6.jpg'    }
%     {'./stimuli/Apfelessen1.jpg'       }
%     {'./stimuli/Apfelessen2.jpg'       }
%     {'./stimuli/Apfelessen3.jpg'       }
%     {'./stimuli/Apfelessen4.jpg'       }
%     {'./stimuli/Apfelessen5.jpg'       }
%     {'./stimuli/Apfelessen6.jpg'       }
%     {'./stimuli/Kuchenessen1.jpg'      }
%     {'./stimuli/Kuchenessen2.jpg'      }
%     {'./stimuli/Kuchenessen3.jpg'      }
%     {'./stimuli/Kuchenessen4.jpg'      }
%     {'./stimuli/Kuchenessen5.jpg'      }
%     {'./stimuli/Kuchenessen6.jpg'      }
%     {'./stimuli/Fensterputzen1.jpg'    }
%     {'./stimuli/Fensterputzen2.jpg'    }
%     {'./stimuli/Fensterputzen3.jpg'    }
%     {'./stimuli/Fensterputzen4.jpg'    }
%     {'./stimuli/Fensterputzen5.jpg'    }
%     {'./stimuli/Fensterputzen6.jpg'    }
%     {'./stimuli/Geschirrabwaschen1.jpg'}
%     {'./stimuli/Geschirrabwaschen2.jpg'}
%     {'./stimuli/Geschirrabwaschen3.jpg'}
%     {'./stimuli/Geschirrabwaschen4.jpg'}
%     {'./stimuli/Geschirrabwaschen5.jpg'}
%     {'./stimuli/Geschirrabwaschen6.jpg'}
%     {'./stimuli/Zähneputzen1.jpg'      }
%     {'./stimuli/Zähneputzen2.jpg'      }
%     {'./stimuli/Zähneputzen3.jpg'      }
%     {'./stimuli/Zähneputzen4.jpg'      }
%     {'./stimuli/Zähneputzen5.jpg'      }
%     {'./stimuli/Zähneputzen6.jpg'      }
%     {'./stimuli/Gesichtwaschen1.jpg'   }
%     {'./stimuli/Gesichtwaschen2.jpg'   }
%     {'./stimuli/Gesichtwaschen3.jpg'   }
%     {'./stimuli/Gesichtwaschen4.jpg'   }
%     {'./stimuli/Gesichtwaschen5.jpg'   }
%     {'./stimuli/Gesichtwaschen6.jpg'   }
