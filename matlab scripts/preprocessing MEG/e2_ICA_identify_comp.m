%% VS + RFT
% PhD project 2

% e. Run ICA with maximum 68 components - reject components

% [c] Katharina Duecker

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions

clear all; close all; clc; beep off
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
%addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA projections
icapth = fullfile(pth,'results','meg', '3 ICA', '1 all subj');
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(icapth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% identify bad components
for s = 15:length(subjfolds)
    load(fullfile(icapth,subjfolds{s}))
    
    cfg = [];
    cfg.layout = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    ft_databrowser(cfg, dataICA)
    shg
    pause
    q = input('continue? y/n','s');

    while ~strcmp(q,'y')
        q = input('continue? y/n','s');
        pause
    end
    
    badcomps = input('bad components? []');
    save(fullfile(icapth,subjfolds{s}), 'badcomps','-append')

    clear dataICA
    close all
end