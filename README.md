# Alpha oscillations support the efficiency of guided visual search by inhibiting both target and distractor features in early visual cortex

The repository contains all scripts associated with a Magnetoencephalography (MEG) - Rapid Invisible Frequency Tagging study aimed at investigating the neural correlates of Guided Search.

authors: Katharina Duecker, Kimron L. Shapiro, Simon Hanslmayr, Jeremy Wolfe, Yali Pan, and Ole Jensen

## MNE scripts

### Maxfilter

Filtering based on maxwell equations and denoising of the sensors is implemented in MNE python.


## Experiment code

experiment/a_exp is the main experiment file;
only works with Propixx Lite projector and has been developed for MATLAB/2017a


## MATLAB scripts

Used for all analyses in the manuscript

### alpha
- a: TFR 4 - 30 Hz for each participant
- b: find Individual Alpha Frequency and sensors of interest
- c: align TFR to IAF
- d: alpha power per condition; performance ~ median split high/low alpha
- e: confirmatory analysis TFR for fast vs slow trials

### behaviour
- a: behaviour per condition
- b: behaviour for alpha high/low
- c: confirmatory analysis: correlation reaction time ~ alpha power (mentioned in manuscript but not explicitly plotted)

### cbrewer

brewer colormap

### coherence
- a: coherence collapsed over trials for topo plots
- b: coherence per condition (priority map)
- c: coherence for fast vs slow (only weak relationship)
- d: coherence for alpha high vs low


### eye movement
- a: ocular artefacts for alpha high vs low
- b: ocular artefacts for fast vs slow trials

### preprocessing MEG
- a: merge .fif; .edf; and .mat (behaviour) files
- b: semi-automatic artefact rejection
- c: find the delay between trigger and diode onset (important to replace diode signal with sine wave)
- d: identify saccades (this is exploratory, there were too many saccades so we controlled for them with analyses mentioned above)
- e: ICA
- f: find tagging sensors of interest -> decide which stimuli to keep
- g: separate data into conditions (this saves one .mat files with MEG data and performance for each condition)

### source: DICS beamformer
- a: align T1 to digitized fiducials
- b: estimate leadfield
- c: DICS for RIFT, alpha pre-search and during search
- d: plot
