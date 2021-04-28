# eegmatlab

Toolboxes required: FieldTrip, EELGAB (plugins ASR, ERPLAB). Data recorded using BioSemi64. MATLAB 2018b.

Preprocessing in EEGLAB using ASR, line noise filtering, ICA, ERPLAB and analysis in FieldTrip using non-parametric cluster based permutation tests. 

Goal is to analyse spontaneous EEG recordings with variable length event durations. By default the variable length durations are sub-divided into 2s non-overlapping windows to facillitate FFT for analysis. Code (with tweaking) supports both between- and within-group designs.
