# Visual Search & Rapid Frequency Tagging

The repository contains all scripts associated with an Magnetoencephalography (MEG) study aimed at investigating the neural correlates of Guided Search.

Please click [here](https://github.com/katduecker/posters/blob/main/duecker_et_al_ICON_2022.pdf) for the most recent poster on the project.

The analyses is still in progress and will be updated irregularly.

## MNE scripts

### Maxfilter

Filtering based on maxwell equations and denoising of the sensors is implemented in MNE python.

## MATLAB scripts

### Experiment code

implemented in Psychtoolbox

### MEG pre-processing

- Alignment of MEG, behavioral and eye movement data
- Semi-automatic artefact rejection
- Independent Component Analysis
- Identification of sensors showing significant RFT response

all implemented in MATLAB using the fieldtrip toolbox

### MEG spectral analysis

- Time-Frequency representation of power (for frequencies < 30 Hz)
- Coherence between MEG sensors of interest and tagging frequency
