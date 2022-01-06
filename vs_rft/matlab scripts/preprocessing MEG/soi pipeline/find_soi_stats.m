%% VS + RFT
% PhD project 2

% use statistical significance as threshold for soi

% [c] Katharina Duecker

%clear all; close all; clc; beep off

function find_soi_stats(s)
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
%pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
rmpath(genpath('/rds/projects/2018/jenseno_entrainment/fieldtrip'))
% addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% merged file path
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb', 'soi');
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

% list subj that have a merge file
d = dir(mergepth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

mkdir(fullfile(cohpth,subjfolds{s}))

% load clean trial structure
load(fullfile(mergepth, 'docu_merge.mat'))
empty_idx = logical(cell2mat(cellfun(@isempty,mergesubj(:,2),'UniformOutput',false)));
mergesubj = mergesubj(~empty_idx,:)

load(fullfile(mergepth, subjfolds{s},'trl_overlap_meg_el_rsp.mat'))

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
   trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*2.5,meginfo.alltrl_bl{p}(:,3)+fs*2,repmat(-2.5*fs,size(meginfo.alltrl_bl{p},1),1)];
   trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;
    
end

% list fif files
d = dir(fullfile(dtpth,subjfolds{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load in maxfiltered data using trial structure

for p = 1:length(f)
    cfg = [];
    cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = {'MEG'};
    % load data and diode separately (to be able to reject ICA)
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);
    % diodes
    cfg.channel = {'MISC004', 'MISC005'};
    diodes_trl{p}= ft_preprocessing(cfg);

end

data = ft_appenddata([],dtprt{:});
diodes = ft_appenddata([], diodes_trl{:});

% weird trials (subject 15 - first trial messed up
time_end = cell2mat(cellfun(@(x) x(end),diodes.time,'UniformOutput',false));

meginfo.rejtrl_all = [find(meginfo.rejtrl_all);find(time_end~=2)'];

save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'elinfo','rspinfo','meginfo','-v7.3')

% fix sampleinfo
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

% select trials to be kept
if ~isempty(meginfo.rejtrl_all)
    cfg = [];
    cfg.trials = ~ismember(1:length(data.sampleinfo),meginfo.rejtrl_all);
    data = ft_selectdata(cfg,data);  
    diodes = ft_selectdata(cfg,diodes);

    meginfo.alltrl_list(meginfo.rejtrl_all,:) = [];
end


clear dtprt diodes_trl

% average grad positions
d = dir(fullfile(megpth,subjfolds{s}, 'meg'));
files = {d.name};
files = files(strncmp(files,'part',4));

grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(megpth,subjfolds{s},'meg',files{fl}))];
end

% mean channel position
mGrad = grad(1);
% average grad structure
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
clear grad

% concatenate cleaned data and diodes
data = ft_appenddata([],data,diodes);
clear diodes diodes_trl dataclean

% relabel diode: 60 Hz vs 67 Hz
load(fullfile(pth,'experiment','trigdef.mat'))


% trials in which yellow (MISC004) was tagged at 60 Hz

% find trigger
misc4_60_trig = cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'126067'),'UniformOutput',false))...
    + cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'216760'),'UniformOutput',false));

% find trial indices with these triggers
misc4_60_trl = ismember(meginfo.alltrl_list(:,1),[trigdef{logical(misc4_60_trig),1}]);
% trials in which cyan (MISC005) was tagged at 60 Hz
misc5_60_trl = ~misc4_60_trl;

% check if numbers make sense
sum(misc4_60_trl)
sum(misc5_60_trl)

% select data
cfg = [];
cfg.trials = misc4_60_trl;
dat_misc4_60 = ft_selectdata(cfg,data);
cfg.trials = misc5_60_trl;
dat_misc5_60 = ft_selectdata(cfg,data);

clear data

% select diodes and MEG
cfg = [];
cfg.channel = 'MEG';
dat_misc4_60_meg = ft_selectdata(cfg,dat_misc4_60);
dat_misc5_60_meg = ft_selectdata(cfg,dat_misc5_60);
cfg.channel = 'MISC004';
dat_misc4_60_m4 = ft_selectdata(cfg,dat_misc4_60);
dat_misc5_60_m4 = ft_selectdata(cfg,dat_misc5_60);
cfg.channel = 'MISC005';
dat_misc4_60_m5 = ft_selectdata(cfg,dat_misc4_60);
dat_misc5_60_m5 = ft_selectdata(cfg,dat_misc5_60);

clear dat_misc4_60 dat_misc5_60
% relabel
dat_misc4_60_m4.label = {'diode 60'};
dat_misc4_60_m5.label = {'diode 67'};

dat_misc5_60_m5.label = {'diode 60'};
dat_misc5_60_m4.label = {'diode 67'};

% append
dat_misc4_60_newlbl = ft_appenddata([],dat_misc4_60_meg,dat_misc4_60_m4,dat_misc4_60_m5);
dat_misc5_60_newlbl = ft_appenddata([],dat_misc5_60_meg,dat_misc5_60_m5,dat_misc5_60_m4);

clear dat_misc4_60_m4 dat_misc4_60_m5 dat_misc4_60_meg dat_misc5_60_m4 dat_misc5_60_m5 dat_misc5_60_meg
% append
data = ft_appenddata([],dat_misc4_60_newlbl,dat_misc5_60_newlbl);

% separate baseline and flicker
cfg = [];
cfg.latency = [-1 0-1/1000];
data_bsl = ft_selectdata(cfg,data);
cfg.latency = [0.5 1.5-1/1000];
data_rft = ft_selectdata(cfg,data);

clear data
% fourier spectra
cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = [60, 67];
%cfg.pad          = 'nextpow2';
fourier          = ft_freqanalysis(cfg,data_rft);
fourier_bl       = ft_freqanalysis(cfg,data_bsl);

clear data*
% statistical test for each occipital sensor
load(fullfile(pth,'matlab scripts','preprocessing MEG','soi pipeline','occi_sens.mat'))
nchancom = length(occi_soi);
stat_mask_60 = zeros(nchancom,1);

for i = 1:nchancom
    cfg            = [];
    cfg.channel    = {occi_soi{i},'diode 60'};
    fourier_tmp    = ft_selectdata(cfg, fourier);
    fourier_bl_tmp = ft_selectdata(cfg, fourier_bl);
    
    cfg                  = [];
    cfg.parameter        = 'fourierspctrm';
    cfg.statistic        = 'ft_statfun_indepsamplesZcoh';  %%%% take fourierspctrm as input, so no time domain information
    cfg.method           = 'montecarlo';
    cfg.frequency        = 60;
    cfg.tail             = 1; %% right sided, grp1 is bigger than grp2
    cfg.alpha            = 0.01;
    cfg.numrandomization = 10000;
    ntrl_1 = size(fourier.fourierspctrm,1);
    ntrl_2 = size(fourier_bl.fourierspctrm,1);
    design = zeros(1, ntrl_1 + ntrl_2);
    design(1,1:ntrl_1) = 1;
    design(1,(ntrl_1 +1):(ntrl_1 + ntrl_2))= 2;
    cfg.design = design;
    cfg.ivar   = 1;
    cfg.design = design;
    stat = ft_freqstatistics(cfg, fourier_tmp, fourier_bl_tmp);
    stat_mask_60(i) = stat.mask;
    stat_p_60(i,1)  = stat.prob;
    clear fourier_tmp fourier_bl_tmp
end

% save sois for all p

%soi_stat = cell(60,1);

%load(fullfile(pth,'matlab scripts', 'preprocessing','soi_stat_allp.mat'))
soi_stat = occi_soi(logical(stat_mask_60));

save(fullfile(cohpth,subjfolds{s},'soi_stat.mat'),'soi_stat')

