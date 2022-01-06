%% Sanity check scripts Visual Search
% explore coherence

clear all; close all; clc; beep off;

% settings
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'data');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts');
figpth = fullfile(mpth,'results','plots','sanity');
pthout = fullfile(mpth,'results','meg','0 sanity');
mkdir(figpth)
addpath(genpath(mtlpth))
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
%addpath('Z:\fieldtrip')
ft_defaults;
subjcode = 'b3ec';
soi = 'MEG2122';

d = dir(megpth);
folds = {d.name};
folds = folds(strncmp(folds,'202',2));
f = cell2mat(cellfun(@(x) strcmp(x(end-3:end),subjcode),folds,'UniformOutput',false));
subpth = fullfile(megpth,folds{f},'meg');
d = dir(subpth);
files = {d.name};
files(1:2) = [];
files = files(cell2mat(cellfun(@(x) strcmp(x(end-6:end),'sss.fif'),files,'UniformOutput',false)));

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
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 52;
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

%% Coherence matlab functions

frqwdth = 2;                  % width passband     
N = 4;                       % length/order FIR filter
foi = 54:80;
[coh60misc4, psdmeg60misc4, psdmisc60misc4, csd60misc4] = kd_coh_hilb_mtlb(data60misc4, foi, 'MISC004', data60misc4.fsample, N,8,2);


save(fullfile(pthout,[subjcode,'_coh_mtlb_freqbands_N',num2str(N),'_sss.mat']),'coh60misc4', 'psdmeg60misc4', 'psdmisc60misc4','csd60misc4')


%% Plot

% PSD
fig = figure;
[~,soiidx] = intersect(data60misc4.label,soi);
subplot(211)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(psdmisc60misc4(1,:,:)))
title('PSD MISC')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
subplot(212)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(psdmeg60misc4(soiidx,:,:)))
title('PSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
print(fig, fullfile(figpth,['psd_meg_diode_mtlb_filtN_',N,'_sss']),'-dpng')

% coherence & CSD% plot coherence and CSD
fig = figure;
[~,soiidx] = intersect(data60misc4.label,soi);
subplot(211)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(csd60misc4(1,:,:)))
title('CSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
subplot(212)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(coh60misc4(soiidx,:,:)))
title('coherence 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
print(fig, fullfile(figpth,['csd_coh_mtlb_filtN_',num2str(N)]),'-dpng')

