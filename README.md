# Action Hierarchies
The data is from Tonghe Zhuang's PhD thesis on subordinate, basic and superordinate level categories of observed actions. We re-do the representational similarity analysis to identify differences in the time-courses for the decoding.

## Analysis steps
- **run01_ShortToLongEpoch.m**: Re-segment EEG data into epochs of -1.2 to 1.8 s for time-frequency (TFR) analysis. Condition codes are added from stimulus protocol files. An ICA is performed on the segmented data and ICs are saved in folder *../data/ica*.
- **run02_ICAToEEG.m**: The script is applied individually on each participant's IC data. Artifactual ICs are excluded and the remaining ICs are back-transformed in to EEG. Results are saved in folder *../data/cleaned*.
- **run03_TLToERP.m**, **run03_TLToTFR.m**, **run03_TLToITC.m**: Cleaned data is transformed to ERP, TFR, or inter-trial coherence (ITC) per condition, i. e., per image. ERP epoch are re-defined to -0.2 to 0.8s length. Files are saved in folders *../data/erp*, *../data/tfr*, and *../data/itc*. The numbers in column 4 of field trialinfo are the original condition codes from the stimulus protocol files. The numbers in column 5 of field trialinfo are the re-sorted condition codes that match with the row and column indices of the model RDMs. Use these for RSA.
- **run04_decodingERP.m**, **run04_decodingTFR.m**: Binary decoding performed on ERP and TFR data. Results are saved in *../data/DAerp* (72 x 72 x 501) and *../data/DAtfr* (72 x 72 x 27 x 90), respectively.  The sorting in the first two dimensions (72 x 72) correspond to that of the model RDMs.
- 
- run05_stats_decodingTF.m finds significant connected bins in the time-and-frequency data. Results are printed at the console and saved in folder ./data/DA/clustermat.mat. The analogous file fpr ERPs run05b_stats_decodingERP.m does nothing?
- run06_plot_decodingERP.m and run06_plot_decodingERP.m plot the mean decoding accuracy over time or time-frequency bins, respectively.
- run07*_RSA*.m does RSA for TF cluster data (unsorted) and ERP data in many combinations of Tonghe / Gregor segments and sorted / unsorted conditions.
- run08_RSAperbin.m is on TF data, but not on the previously defined clusters, but on each time-frequency bin.
- run09_ClusterCorrectionPerBin.m does post-hoc cluster correction in results from run08_RSAperbin.m

## Open questions
- Thonge liest die Bedingungscodes im script ActiHierEEG/Progs/PairwiseDecoding.m ein, aus Dateien ActiHierEEG/Data/log/SUBxy_stimuli_preprocess.mat. Diese Dateien liegen auf der NAS nicht vor.
- Ich habe auch nach dem Umsortieren in DAerp2 weiterhin andere Ergebnisse als Tonghe.
