%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e2. compare TFRs for fast slow w/ cluster-based test (Fig. 4 d)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs
% - s: subject index
% - c_idx: condition index
% - data_trim (bool): discard trials w/ RT +- 3*std?
% - split_ta_tp (bool): split into target absent/present?

% Output
% - soi_grad: gradiometers with high alpha power
% iaf_grad: identified IAF in gradiometers

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 23/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

clear all; close all; clc; beep off;
condi = {{'ti','16t'}, {'ni','16t'},{'ti','32t'},{'ni','32t'}};

%% settings
pth = 'Z:\Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha','RT');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi/');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','alpha high low');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphafigpth = fullfile(pth,'results','meg','6 Alpha','fig');

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

addpath('Z:\fieldtrip')
ft_defaults;

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));
% load template TFR
% load in example subject to get ft structure
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
fs = 1000;
% list subjects
d = dir(cohpth);
d = {d.name};
subj = d(strncmp(d,'202',3));
d = dir(fullfile(maxfpth,subj{1}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(maxfpth,subj{1},f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = 'MEGGRAD';
% load in data for this part
data = ft_preprocessing(cfg);

winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEG';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'no';
TFR_cond = ft_freqanalysis(cfg,data);

cfg = [];
cfg.method = 'sum';
TFR_cond= ft_combineplanar(cfg,TFR_cond);

cfg = [];
cfg.frequency = [TFR_cond.freq(1), TFR_cond.freq(12)];
TFR_cond = ft_selectdata(cfg,TFR_cond);

pow_fast = cell(1,length(subj));
pow_slow = cell(1,length(subj));
pow_contr = cell(1,length(subj));

% target absent
condi = {{'ti','16ta'}, {'ni','16ta'},{'ti','32ta'},{'ni','32ta'}};

fig = figure('Position',[0 0 1920/1.5 1080]);

for c = 1:length(condi)

    for s = 1:length(subj)

        load(fullfile(outpth,subj{s},[strjoin(condi{c},'_'),'_RTsplit.mat']))
        %TFR_cond.freq = freqvec;

        pow_fast{s} = TFR_cond;
        pow_slow{s} = TFR_cond;

        pow_fast{s}.powspctrm(1,:,:) = iaf_fast;
        pow_slow{s}.powspctrm(1,:,:) = iaf_slow;



        cfg = [];
        cfg.latency = [-1 0.5];
        pow_fast{s} = ft_selectdata(cfg,pow_fast{s});
        pow_slow{s} = ft_selectdata(cfg,pow_slow{s});


        pow_contr{s} = pow_fast{s};

        pow_contr{s}.powspctrm = pow_fast{s}.powspctrm./pow_slow{s}.powspctrm-1;

        cfg = [];
        cfg.channel = TFR_cond.label{1};
        pow_contr{s} = ft_selectdata(cfg,pow_contr{s});
    end

    cfg = [];
    cfg.channel          = TFR_cond.label{1};
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;
    cfg.correcttail  = 'alpha';
    % specifies with which sensors other sensors can form clusters

    design = ones(2,2*length(subj));
    design(1,:) = [1:length(subj), 1:length(subj)];
    design(2,length(subj)+1:2*length(subj)) = 2;

    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;

    [stat] = ft_freqstatistics(cfg, pow_fast{:}, pow_slow{:});

    GA_contr = ft_freqgrandaverage([],pow_contr{:});

    GA_contr.mask = stat.mask;

    stat.freq = freqvec;
    GA_contr.freq = freqvec;

    cfg = [];
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
   %cfg.parameter = 'stat';
    cfg.colormap = cm;
    cfg.zlim = [-0.2 0.2];
    cfg.title = strjoin(condi{c},' ');
    cfg.channel = stat.label{1};
    subplot(2,2,c)
    ft_singleplotTFR(cfg,GA_contr);
    xlabel('time (s)')
    xticks([-1:0.5:0.5])
    ylim([round(freqvec(1)), round(freqvec(end))])
    yticks(round(freqvec(1)):4:round(freqvec(end)))
    cb = colorbar;
    cb.FontName = 'Arial';
    cb.Ticks = -0.2:0.2:0.2;
    
end

print(fig,fullfile(alphafigpth,'alpha_TFR_fast_vs_slow_ta'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'alpha_TFR_fast_vs_slow_ta'),'-dpng','-r0')


% target present
condi = {{'ti','16tp'}, {'ni','16tp'},{'ti','32tp'},{'ni','32tp'}};
fig = figure('Position',[0 0 1920/1.5 1080]);

for c = 1:length(condi)

    for s = 1:length(subj)

        load(fullfile(outpth,subj{s},[strjoin(condi{c},'_'),'_RTsplit.mat']))
        %TFR_cond.freq = freqvec;

        pow_fast{s} = TFR_cond;
        pow_slow{s} = TFR_cond;

        pow_fast{s}.powspctrm(1,:,:) = iaf_fast;
        pow_slow{s}.powspctrm(1,:,:) = iaf_slow;



        cfg = [];
        cfg.latency = [-1 0.5];
        pow_fast{s} = ft_selectdata(cfg,pow_fast{s});
        pow_slow{s} = ft_selectdata(cfg,pow_slow{s});


        pow_contr{s} = pow_fast{s};

        pow_contr{s}.powspctrm = pow_fast{s}.powspctrm./pow_slow{s}.powspctrm-1;
    end

    cfg = [];
    cfg.channel          = TFR_cond.label{1};
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;
        cfg.correcttail  = 'alpha';

    % specifies with which sensors other sensors can form clusters

    design = ones(2,2*length(subj));
    design(1,:) = [1:length(subj), 1:length(subj)];
    design(2,length(subj)+1:2*length(subj)) = 2;

    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;

    [stat] = ft_freqstatistics(cfg, pow_fast{:}, pow_slow{:});

    GA_contr = ft_freqgrandaverage([],pow_contr{:});

    GA_contr.mask = stat.mask;

    stat.freq = freqvec;
    GA_contr.freq = freqvec;

    cfg = [];
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
   %cfg.parameter = 'stat';
    cfg.colormap = cm;
    cfg.zlim = [-0.2 0.2];
    cfg.title = strjoin(condi{c},' ');
    cfg.channel = stat.label{1};
    subplot(2,2,c)
    ft_singleplotTFR(cfg,GA_contr);
    xlabel('time (s)')
    xticks([-1:0.5:0.5])
    cb = colorbar;
    cb.FontName = 'Arial';
    cb.Ticks = -0.2:0.2:0.2;
    ylim([round(freqvec(1)), round(freqvec(end))])
    yticks(round(freqvec(1)):4:round(freqvec(end)))
    
end

print(fig,fullfile(alphafigpth,'alpha_TFR_fast_vs_slow_tp'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'alpha_TFR_fast_vs_slow_tp'),'-dpng','-r0')

