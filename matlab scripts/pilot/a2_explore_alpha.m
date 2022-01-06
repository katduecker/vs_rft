%% TFR in low frequency bands
clear all; close all; clc; beep off;

% settings
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'data');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts');
figpth = fullfile(mpth,'results','plots','sanity');
mkdir(figpth)
addpath(genpath(mtlpth))
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
%addpath('Z:\fieldtrip')
ft_defaults;
subjcode = 'b583';

d = dir(megpth);
folds = {d.name};
folds = folds(strncmp(folds,'202',2));
f = cell2mat(cellfun(@(x) strcmp(x(end-3:end),subjcode),folds,'UniformOutput',false));
subpth = fullfile(megpth,folds{f},'meg');
d = dir(subpth);
files = {d.name};
files = files(strncmp(files,'part',4));

clear d f folds

% average grad positions

% load in grad structures
grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(subpth,files{fl}))];
end
% average
mGrad = grad(1);
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
clear grad

%% Load in data

% load trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(subpth,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    % target or distractor tagging?
    cfg.detrend = 'yes';
    %cfg.hpfilter = 'yes';
    %cfg.hpfreq = 52;
    % load in all trials
    cfg.trl = kd_trlfun_alltrl(cfg); 
    block{fl} = ft_preprocessing(cfg);
end

data = ft_appenddata([],block{:});
fs = block{fl}.fsample;                 % sampling rate
clear block

% pad data to 8 sec
% zero pad before averaging
t = -2.5:1/fs:5.5-1/fs;
datapad = kd_datapad(data,fs,8,t);

% phase-locked alpha
cfg = [];
ERF = ft_timelockanalysis(cfg,datapad);

%% TFR alpha band

winl = 0.5;
% frequency analysis
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
% cfg.pad = 8;
% cfg.padtype = 'mean';
cfg.foi = 4:1:30;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
cfg.toi = -1.25:0.05:1.25;
cfg.keeptrials = 'no';
TFRphlock = ft_freqanalysis(cfg,ERF);  
TFRspont = ft_freqanalysis(cfg,datapad);  
TFRphlock.grad = mGrad;
TFRspont.grad = mGrad;

% combine
cfg = [];
cfg.method = 'sum';
TFRphlock = ft_combineplanar(cfg,TFRphlock);
TFRspont = ft_combineplanar(cfg,TFRspont);

% TFRs
fig = figure;
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.colorbar = 'yes';
cfg.channel = {'MEG2012+2013'};
subplot(211)
ft_singleplotTFR(cfg,TFRphlock);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('all trials soi MEG2012+2013; phase locked alpha')
subplot(212)
ft_singleplotTFR(cfg,TFRspont);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('all trials soi MEG2012+2013; non-phase locked alpha')
print(fig,fullfile(figpth,[subjcode,'_alpha_tfr_',num2str(winl*1000)]),'-dpng')

% average -> power spectra
cfg = [];
cfg.latency = [-1.5 -winl/2];
cfg.avgovertime = 'yes';
AVGprelock = ft_selectdata(cfg,TFRphlock);
AVGprenlock = ft_selectdata(cfg,TFRspont);
cfg.latency = [winl/2 1.5];
cfg.avgovertime = 'yes';
AVGperilock = ft_selectdata(cfg,TFRphlock);
AVGperinlock = ft_selectdata(cfg,TFRspont);

soiidx = find(strcmp(AVGprelock.label,'MEG2012+2013'));
fig = figure;
subplot(221)
plot(AVGprelock.freq,AVGprelock.powspctrm(soiidx,:));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('pre-stim alpha, ph-lock')
subplot(222)
plot(AVGperilock.freq,AVGperilock.powspctrm(soiidx,:));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('peri alpha, ph-lock')
subplot(223)
plot(AVGprenlock.freq,AVGprenlock.powspctrm(soiidx,:));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('pre-stim alpha, non-ph-lock')
subplot(224)
plot(AVGperinlock.freq,AVGperinlock.powspctrm(soiidx,:));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('peri alpha, non-ph-lock')
print(fig,fullfile(figpth,[subjcode,'_alpha_spectra_',num2str(winl*1000)]),'-dpng')
