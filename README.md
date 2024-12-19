# Action Hierarchies
The data is from Tonghe Zhuang's PhD thesis on subordinate, basic and superordinate level categories of observed actions. We re-do the representational similarity analysis to identify differences in the time-courses for the decoding.

## Analysis steps
- run01_ShortToLongEpoch.m Re-segment EEG data into epochs of -1.2 to 1.8 s for time-frequency (TFR) analysis. Condition codes are added from stimulus protocol files. An ICA is performed on the segmented data and ICs are saved in folder ./data/ica.
- run02_ICAToEEG.m This script is applied individually on each participant's IC data. Artifactual ICs are excluded and the remaoning ICs are back-transformed in to EEG. Results are saved in folder ./data/cleaned-
- run03_TLToTFR.m, run03b_TLToITC.m, run03c_TLToERP.m Cleaned data is transformed to the mean TFR, inter-trial coherence (ITC) and ERP per condition, i. e., per image. Files are saved in folders ./data/tfr, ./data/itc and ./data/erp. Mind that the condition codes do NOT correspond to row and column indices of the model RDMs and need to be re-sorted before the RSA (see below).
- NEW: run03_TLToERP.m has re-sorted conditon codes in column 5 of field trialinfo
- run04_decoding.m etc. Binary decoding is performed on TFR, ICA and ERP data. Results are saved in ./data/DA and ./data/DAtfr2 for TFR; and ./data/DAerp, ./data/DAerp2 and ./data/DAerp3 for ERPs. In ./data/DA and ./data/DAerp, condition codes 1:72 correspond to codes 1:72 in the stimulus protocol files. Folder ./data/DAerp2 contains results from Tonghe's original EEG segments and condition codes are re-ordered to match row and column order of the model RDMs. ./data/DA3erp is same as ./data/DAtfr2?
- run05_stats_decodingTF.m finds significant connected bins in the time-and-frequency data. Results are printed at the console and saved in folder ./data/DA/clustermat.mat. The analogous file fpr ERPs run05b_stats_decodingERP.m does nothing?
- run06_plot_decodingERP.m and run06_plot_decodingERP.m plot the mean decoding accuracy over time or time-frequency bins, respectively.
- run07*_RSA*.m does RSA for TF cluster data (unsorted) and ERP data in many combinations of Tonghe / Gregor segments and sorted / unsorted conditions.
- run08_RSAperbin.m is on TF data, but not on the previously defined clusters, but on each time-frequency bin.
- run09_ClusterCorrectionPerBin.m does post-hoc cluster correction in results from run08_RSAperbin.m

## Open questions
- Thonge liest die Bedingungscodes im script ActiHierEEG/Progs/PairwiseDecoding.m ein, aus Dateien ActiHierEEG/Data/log/SUBxy_stimuli_preprocess.mat. Diese Dateien liegen auf der NAS nicht vor.
- Ich habe auch nach dem Umsortieren in DAerp2 weiterhin andere Ergebnisse als Tonghe.
