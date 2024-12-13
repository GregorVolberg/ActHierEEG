%% convert the order of stimuli as the same with fMRI data
% from Tonghe
function X_sorted=EEG_reshape_to_fMRI(X)

orderOfActions=[49,50,51,52,53,54,13,14,15,16,17,18,37,38,39,40,41,42,55,56,...
    57,58,59,60,7,8,9,10,11,12,61,62,63,64,65,66,1,2,3,4,5,6,43,44,45,46,47,48,...
    19,20,21,22,23,24,25,26,27,28,29,30,67,68,69,70,71,72,31,32,33,34,35,36];

X_sorted=X(orderOfActions, orderOfActions);