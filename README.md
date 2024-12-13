# Action Hierarchies
The data is from Tonghe Zhuang's PhD thesis on subordinate, basic and superordinate level categories of observed actions. We re-do the representational similarity analysis to identify differences in the time-courses for the decoding.

## Analysis steps
- run01_ShortToLongEpoch.m Re-segment EEG data into epochs of -1.2 to 1.8 s for time-frequency (TFR) analysis. Condition codes are added from stimulus protocol files. An ICA is performed on the segmented data and ICs are saved in folder ./data/ica.
- run02_ICAToEEG.m This script is applied individually on each participant's IC data. Artifactual ICs are excluded and the remaoning ICs are back-transformed in to EEG. Results are saved in folder ./data/cleaned-
- run03_TLToTFR.m, run03b_TLToITC.m, run03c_TLToERP.m Cleaned data is transformed to the mean TFR, inter-trial coherence (ITC) and ERP per condition, i. e., per image. Files are saved in folders ./data/tfr, ./data/itc and ./data/erp. Mind that the condition codes do NOT correspond to row and column indices of the model RDMs and need to be re-sorted before the RSA (see below).
- run04_decoding.m etc. Binary decoding is performed on TFR, ICA and ERP data. Results are saved in ./data/DA and ./data/DAtfr2 for TFR; and ./data/DAerp, ./data/DAerp2 and ./data/DAerp3 for ERPs. [What is what?]
- run05_stats_decoding.m Does statistics on decoding accuracy.
