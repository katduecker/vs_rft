%% VS + RFT
% PhD project 2

% artefact suppression using ICA
% [c] Katharina Duecker

clear all; close all; clc; beep off
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% 1. read in data in trials +- 1 second
for s = 1:length(subjfolds)
    % load in trial structure
    load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'alltrl_bl')
    
    % trial structure to load in trl
    for p = 1:length(alltrl_bl)
        trlstruct{p} = [alltrl_bl{p}(:,2)-fs,alltrl_bl{p}(:,3)+fs,zeros(length(alltrl_bl{p}),1)];
    end
    
    % load in maxfiltered data using trial structure
    
    % list fif files
    d = dir(fullfile(dtpth,subjfolds{s}));
    f = {d.name};
    % find fif file
    idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    f = f(idxx);
    
    cfg = [];
    for p = 1:length(f)
        cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
        cfg.detrend = 'yes';
        cfg.trl = trlstruct{p};
        cfg.channel = 'MEG'
        % load in data for this part
        dtprt{p} = ft_preprocessing(cfg);
    end
    
    data = ft_appenddata([],dtprt{:});
    
    % fix sample info for downsampling
    
    sinfo = zeros(length(data.trial),1);
    % first trial samples: 0 to length of data in samples
    sinfo(1,2) = size(data.trial{1},2);
    for t = 2:length(data.trial)
        % next trial: one sample after end of previous
        sinfo(t,1) = sinfo(t-1,2) + 1;
        % start sample + length of trial - 1
        sinfo(t,2) = sinfo(t,1) + size(data.trial{t},2);
    end
    
    data.sampleinfo = sinfo;
    
    % resample to 250 Hz (1/4 of sampling rate)
    
    cfg = [];
    cfg.resamplefs = 250;
    dataresamp = ft_resampledata(cfg,data);
    
    % 2. ICA
    cfg = [];
    cfg.method = 'runica';
    dataICA = ft_componentanalysis(cfg,dataresamp);
    
    %mkdir(fullfile(pth,'results','meg','3 ICA'))
    save(fullfile(pth,'results','meg','3 ICA',[subjfolds{s},'_ica.mat']),'dataICA','-v7.3')
    clear data* dtprt
end