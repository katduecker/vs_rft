%% VS + RFT
% PhD project 2

% Pre-processing of usable participants

% 1. Identify eye movement trials & store
% 2. read in trials based on identified trial structure (merged between
% MEG, EDF, RT)
% 3a. semi-automatic artefact rejection
% 3b. overlap trials eye movement & artefact rejection?
% 3c. ICA
% 4a. Split data into eye movement and no movement 
% 4b. Coherence separately for movement and no movement

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path
dtpth = fullfile(pth, 'data');                                         % data path
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

s = 1;                                                                 % subject id

%% Identify eye movement trials & store
a1_fun_identify_saccades(s, pth, dtpth)

%% Semi-automatic artefact rejection
fs = 1000;
a2_fun_artef_rej(s,pth,fs);

%% ICA
