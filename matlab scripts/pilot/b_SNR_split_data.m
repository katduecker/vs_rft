%% Sanity check scripts Visual Search
% SNR when data are split according to Target/Distractor tagging

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
subjcode = 'b3ec';

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
    % separate trials based on target tagging
    cfg.tag = 't';
    cfg.detrend = 'yes';
    %cfg.demean = 'yes';
    %cfg.baselinewindow = [-1.25 -0.25];
%     cfg.hpfilter = 'yes';
%     cfg.hpfreq = 52;
    % load in all trials
    [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg); 
    cfg.trl = trl60_misc4;
    block60misc4{fl} = ft_preprocessing(cfg);
    cfg.trl = trl60_misc5;   
    block60misc5{fl} = ft_preprocessing(cfg);
    cfg.trl = trl67_misc4;
    block67misc4{fl} = ft_preprocessing(cfg);
    cfg.trl = trl67_misc5;
    block67misc5{fl} = ft_preprocessing(cfg);
end

data60misc4 = ft_appenddata([],block60misc4{:});
data60misc4.fsample = block60misc4{1}.fsample;

data60misc5 = ft_appenddata([],block60misc5{:});
data60misc5.fsample = block60misc4{1}.fsample;

data67misc4 = ft_appenddata([],block67misc4{:});
data67misc4.fsample = block60misc4{1}.fsample;

data67misc5 = ft_appenddata([],block67misc5{:});
data67misc5.fsample = block60misc4{1}.fsample;

clear block*

% zero pad before averaging
fs = data60misc4.fsample;
t = -2.5:1/fs:5.5-1/fs;
data60misc4 = kd_datapad(data60misc4,fs,8,t);
data60misc5 = kd_datapad(data60misc5,fs,8,t);
data67misc4 = kd_datapad(data67misc4,fs,8,t);
data67misc5 = kd_datapad(data67misc5,fs,8,t);


% Compute ERF
ERF60m4 = ft_timelockanalysis([],data60misc4);
ERF60m5 = ft_timelockanalysis([],data60misc5);
ERF67m4 = ft_timelockanalysis([],data67misc4);
ERF67m5 = ft_timelockanalysis([],data67misc5);

% % baseline correct
% cfg = [];
% cfg.baseline = [-1.25 -0.25];
% ERF60m4 = ft_timelockbaseline(cfg,ERF60m4);

misc4 = strcmp(ERF60m4.label,'MISC004');
misc5= strcmp(ERF60m5.label,'MISC005');

subplot(421)
plot(ERF60m4.time,ERF60m4.avg(misc4,:))
title('60 Hz')
xlim([-0.25 0.75])
subplot(422)
plot(ERF60m4.time,ERF60m4.avg(misc5,:))
title('67 Hz')
xlim([-.25 0.75])

subplot(423)
plot(ERF60m5.time,ERF60m5.avg(misc4,:))
title('67 Hz')
xlim([-.25 0.75])
subplot(424)
plot(ERF60m5.time,ERF60m5.avg(misc5,:))
title('60 Hz')
xlim([-.25 0.75])

subplot(425)
plot(ERF67m4.time,ERF67m4.avg(misc4,:))
title('67 Hz')
xlim([-.25 0.75])
subplot(426)
plot(ERF67m4.time,ERF67m4.avg(misc5,:))
title('60 Hz')
xlim([-.25 0.75])

subplot(427)
plot(ERF67m5.time,ERF67m5.avg(misc4,:))
title('60 Hz')
xlim([-.25 0.75])
subplot(428)
plot(ERF67m5.time,ERF67m5.avg(misc5,:))
title('67 Hz')
xlim([-.25 0.75])
print(fig,fullfile(figpth,'diode phase-locked'),'-dpng')

%% TFRs

% contrast conditions
winl = 0.5;
% frequency analysis
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
% cfg.pad = 8;
% cfg.padtype = 'mean';
cfg.foi = 4:1:80;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
cfg.toi = -1.25:0.05:1.25;
cfg.keeptrials = 'no';
TFR60misc4 = ft_freqanalysis(cfg,ERF60m4);                  % evoked/phase-locked
TFR67misc4 = ft_freqanalysis(cfg,ERF67m4);                  % evoked/phase-locked
TFR60misc5 = ft_freqanalysis(cfg,ERF60m5);                  % evoked/phase-locked
TFR67misc5 = ft_freqanalysis(cfg,ERF67m5);                  % evoked/phase-locked
TFR60misc4.grad = mGrad;
TFR60misc5.grad = mGrad;
TFR67misc4.grad = mGrad;
TFR67misc5.grad = mGrad;

cfg = [];
cfg.method = 'sum';
TFR60misc4 = ft_combineplanar(cfg,TFR60misc4);
TFR67misc4 = ft_combineplanar(cfg,TFR67misc4);
TFR60misc5 = ft_combineplanar(cfg,TFR60misc5);
TFR67misc5 = ft_combineplanar(cfg,TFR67misc5);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [58 80];
cfg.baseline = [-1.5 -0.5];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
%cfg.zlim = [0 8];
%cfg.xlim = [0.5 4];
%ft_multiplotTFR(cfg,TFR60misc4)

cfg.channel = {'MEG2012+2013'};
fig = figure;
subplot(211)
ft_singleplotTFR(cfg,TFR60misc4);
ylabel('frequency (Hz)')
xlabel('time (s)')
title(['misc4 tagged at 60 Hz ',num2str(winl),'s win; averaged TFR'])

% average TFR
cfg = [];
cfg.baseline = [-1.5 -0.5];
cfg.baselinetype = 'relchange';
TFR60m4bsl = ft_freqbaseline(cfg,TFR60misc4);

cfg = [];
cfg.latency = [0.25 1];
cfg.avgovertime = 'yes';
AVG60m4 = ft_selectdata(cfg,TFR60m4bsl);

soiidx = find(strcmp(AVG60m4.label,'MEG2012+2013'));
subplot(212)
plot(AVG60m4.powspctrm(soiidx,:))
xlim([40 80])
ylabel('frequency (Hz)')
xlabel('time (s)')
title(['misc4 tagged at 60 Hz ',num2str(winl),'s win; averaged TFR'])

print(fig,fullfile(figpth,['TFR_winl_',num2str(winl)]),'-dpng')

%% Realtive power change

soi = {'MEG2032','MEG2033', 'MEG2042'};
soicmb = {'MEG2032+2033', 'MEG2042+2043'};
bsltoi = [-1.2 -0.2];
stimtoi = [0.25 1.25];
foi = [54 80];

% evoked spectra

[~,~,bsl, stim60m4] = kd_fft_bsl_stim(ERF60m4,soi,'MISC004','MISC005',bsltoi,stimtoi,'no',foi,mGrad);

% calculate relative change
RC60m4 = stim60m4;
RC60m4.powspctrm = stim60m4.powspctrm./bsl.powspctrm - 1;

% evoked spectra
[~,~,bsl, stim60m5] = kd_fft_bsl_stim(ERF60m5,soi,'MISC004','MISC005',bsltoi,stimtoi,'no',foi,mGrad);

% calculate relative change
RC60m5 = stim60m5;
RC60m5.powspctrm = stim60m5.powspctrm./bsl.powspctrm - 1;

% evoked spectra
[~,~,bsl, stim67m4] = kd_fft_bsl_stim(ERF67m4,soi,'MISC004','MISC005',bsltoi,stimtoi,'no',foi,mGrad);

% calculate relative change
RC67m4 = stim67m4;
RC67m4.powspctrm = stim67m4.powspctrm./bsl.powspctrm - 1;

% evoked spectra
[~,~,bsl, stim67m5] = kd_fft_bsl_stim(ERF67m5,soi,'MISC004','MISC005',bsltoi,stimtoi,'no',foi,mGrad);

% calculate relative change
RC67m5 = stim67m5;
RC67m5.powspctrm = stim67m5.powspctrm./bsl.powspctrm - 1;

% Tag SOI
fig = figure;
cfg = [];
cfg.channel = soicmb;
subplot(221)
ft_singleplotER(cfg,RC60m4)
title('T misc4 60 Hz')
subplot(222)
ft_singleplotER(cfg,RC60m5)
title('T misc5 60 Hz')
subplot(223)
ft_singleplotER(cfg,RC67m4)
title('T misc4 67 Hz')
subplot(224)
ft_singleplotER(cfg,RC67m5)
title('T misc5 67 Hz')
print(fig,fullfile(figpth,'RC_T_tag_soi_misc'),'-dpng')



fig = figure;
cfg = [];
cfg.channel = 'MISC004';
subplot(221)
ft_singleplotER(cfg,RC60m4)
title('T misc4 60 Hz')
subplot(222)
cfg.channel = 'MISC005';
ft_singleplotER(cfg,RC60m5)
title('T misc5 60 Hz')
subplot(223)
cfg.channel = 'MISC004';
ft_singleplotER(cfg,RC67m4)
title('T misc4 67 Hz')
subplot(224)
cfg.channel = 'MISC005';
ft_singleplotER(cfg,RC67m5)
title('T misc5 67 Hz')
print(fig,fullfile(figpth,'RC_tag_diode'),'-dpng')





%% Separate based on T frequency
% 60 Target vs 60 Distractor
dataT60 = ft_appenddata([],data60misc4,data60misc5);
dataT67 = ft_appenddata([],data67misc4,data67misc5);

% zero pad before averaging
t = -2.5:1/fs:5.5-1/fs;
fs = data60misc4.fsample;
dataT60 = kd_datapad(dataT60,fs,8,t);
dataT67 = kd_datapad(dataT67,fs,8,t);


% % Compute ERF
% ERFT60 = ft_timelockanalysis([],dataT60);
% ERFT67 = ft_timelockanalysis([],dataT67);

% FFT
% evoked
[~,~,bsl, stimT60ev] = kd_fft_bsl_stim(dataT60,soi,'MISC004','MISC005',bsltoi,stimtoi,'yes',foi,mGrad);
RCT60ev = stimT60ev;
RCT60ev.powspctrm = stimT60ev.powspctrm./bsl.powspctrm - 1;
[~,~,bsl, stimT67ev] = kd_fft_bsl_stim(dataT67,soi,'MISC004','MISC005',bsltoi,stimtoi,'yes',foi,mGrad);
RCT67ev = stimT67ev;
RCT67ev.powspctrm = stimT67ev.powspctrm./bsl.powspctrm - 1;

% induced
[~,~,bsl, stimT60ind] = kd_fft_bsl_stim(dataT60,soi,'MISC004','MISC005',bsltoi,stimtoi,'no',foi,mGrad);
RCT60ind = stimT60ind;
RCT60ind.powspctrm = stimT60ind.powspctrm./bsl.powspctrm - 1;
[~,~,bsl, stimT67ind] = kd_fft_bsl_stim(dataT67,soi,'MISC004','MISC005',bsltoi,stimtoi,'no',foi,mGrad);
RCT67ind = stimT67ind;
RCT67ind.powspctrm = stimT67ind.powspctrm./bsl.powspctrm - 1;

% plot power change
fig = figure;
subplot(221)
cfg = [];
cfg.channel = soicmb;
ft_singleplotER(cfg,stimT60ev);
title('T tagged at 60 Hz evoked')
ylabel('power (T/m)²')
subplot(222)
ft_singleplotER(cfg,RCT60ev);
title('T tagged at 60 Hz evoked')
ylabel('rel chan')
subplot(223)
ft_singleplotER(cfg,stimT67ev);
title('T tagged at 67 Hz evoked')
ylabel('power (T/m)²')
subplot(224)
ft_singleplotER(cfg,RCT67ev);
title('T tagged at 67 Hz evoked')
ylabel('rel chan')
print(fig,fullfile(figpth,'FFT_T_tag_ev'),'-dpng')

fig = figure;
subplot(221)
cfg = [];
cfg.channel = soicmb{1};
ft_singleplotER(cfg,stimT60ev);
title('T tagged at 60 Hz evoked')
ylabel('power (T/m)²')
subplot(222)
ft_singleplotER(cfg,RCT60ev);
title('T tagged at 60 Hz evoked')
ylabel('rel chan')
subplot(223)
ft_singleplotER(cfg,stimT67ev);
title('T tagged at 67 Hz evoked')
ylabel('power (T/m)²')
subplot(224)
ft_singleplotER(cfg,RCT67ev);
title('T tagged at 67 Hz evoked')
ylabel('rel chan')
print(fig,fullfile(figpth,['FFT_T_tag_ev_',soicmb{1}]),'-dpng')

fig = figure;
subplot(221)
cfg = [];
cfg.channel = soicmb{2};
ft_singleplotER(cfg,stimT60ev);
title('T tagged at 60 Hz evoked')
ylabel('power (T/m)²')
subplot(222)
ft_singleplotER(cfg,RCT60ev);
title('T tagged at 60 Hz evoked')
ylabel('rel chan')
subplot(223)
ft_singleplotER(cfg,stimT67ev);
title('T tagged at 67 Hz evoked')
ylabel('power (T/m)²')
subplot(224)
ft_singleplotER(cfg,RCT67ev);
title('T tagged at 67 Hz evoked')
ylabel('rel chan')
print(fig,fullfile(figpth,['FFT_T_tag_ev_',soicmb{2}]),'-dpng')

% plot power change
fig = figure;
subplot(221)
cfg = [];
cfg.channel = soicmb;
ft_singleplotER(cfg,stimT60ev);
title('T tagged at 60 Hz evoked')
subplot(222)
cfg = [];
cfg.channel = soicmb;
ft_singleplotER(cfg,RCT60ev);
title('T tagged at 60 Hz evoked')
ylabel('power change relative to bsl')
subplot(223)
cfg = [];
cfg.channel = soicmb;
ft_singleplotER(cfg,stimT67ev);
title('T tagged at 67 Hz evoked')
subplot(224)
cfg = [];
cfg.channel = soicmb;
ft_singleplotER(cfg,RCT67ev);
title('T tagged at 67 Hz evoked')
ylabel('power change relative to bsl')
print(fig,fullfile(figpth,'FFT_T_tag_ev'),'-dpng')


% Contrast conditions
Tenh60 = stimT60ev;
Tenh60.powsptrm = stimT60ev.powspctrm./stimT67ev.powspctrm - 1;

Tenh67 = stimT67ev;
Tenh67.powspctrm = stimT67ev.powspctrm./stimT60ev.powspctrm - 1;

fig = figure;
cfg.channel = soicmb;
ft_singleplotER(cfg,Tenh60);
print(fig,fullfile(figpth,'Contrast conditions T vs D tagging'),'-dpng')
