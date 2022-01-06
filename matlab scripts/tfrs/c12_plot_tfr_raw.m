%% VS + RFT
% PhD project 2

% plots TFR of power on raw data

% [c] Katharina Duecker

clear all; close all; clc; beep off;
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
tfrpth = fullfile(pth,'results','meg', '4 TFR power','raw');

% figure paths
figpth = fullfile(pth,'results','meg', 'figures', 'TFR power');
mkdir(figpth)
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(megpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

s = 8;

load('soi_tfr_subj.mat')

load(fullfile(tfrpth,[subjfolds{s},'_TFRalltrl.mat']))
%% Plots

% multiplot
% MAG
figure('Name', 'MAG')
cfg = [];
cfg.xlim = [-1.5 1];
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
cfg.channel = 'MEGMAG';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'yes';
ft_multiplotTFR(cfg,TFRall)
%soimag{s} = {'MEG1911', 'MEG2031', 'MEG2041'};

% GRAD
cfg = [];
cfg.method = 'sum';
TFRall = ft_combineplanar(cfg,TFRall);

figure('Name', 'GRAD')
cfg = [];
cfg.xlim = [-1.5 1];
cfg.ylim = [56 72];
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
cfg.channel = 'MEGGRAD';
cfg.layout = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg,TFRall)
%soicmb{s} = {'MEG1932+1933', 'MEG2032+2033', 'MEG2042+2043', 'MEG2112+2113', 'MEG2342+2343'};

% save soi
%save('soi_tfr_subj.mat','soimag','soicmb')
fig = figure('Name', 'TFR avg SOI');
cfg = [];
cfg.xlim = [-1.25 1];
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
cfg.channel = soimag{s};
% cfg.xlim = [-1 2];
% cfg.ylim = [56 74];
cfg.title = 'magnetometers OI';
subplot(211)
ft_singleplotTFR(cfg,TFRall);
cb = colorbar;
cb.Label.String = 'relative power change';
cfg.title = 'cmb grads OI';
line(TFRall.time,repmat(67,1,length(TFRall.time)),'Color','k','LineStyle', '-.')
line(TFRall.time,repmat(60,1,length(TFRall.time)),'Color','k','LineStyle', '-.')
xlim([-1.25 1])
xlabel('time (s)')
ylabel('frequency (Hz)')
subplot(212)
cfg.channel = soicmb{s};
ft_singleplotTFR(cfg,TFRall);
cb = colorbar;
cb.Label.String = 'relative power change';
line(TFRall.time,repmat(67,1,length(TFRall.time)),'Color','k','LineStyle', '-.')
line(TFRall.time,repmat(60,1,length(TFRall.time)),'Color','k','LineStyle', '-.')
xlim([-1.25 1])
xlabel('time (s)')
ylabel('frequency (Hz)')
print(fig,fullfile(figpth,subjfolds{s},'tfr_mags_cmb_raw'),'-dpng')
close all

% apply baseline and average over time
cfg = [];
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
TFRall = ft_freqbaseline(cfg,TFRall);

cfg = [];
cfg.latency = [.5 1];
cfg.avgovertime = 'yes';
AVGPOW = ft_selectdata(cfg, TFRall);

% plot
fig = figure('Name', 'avg TFR avg SOI');
subplot(211)
plot(AVGPOW.freq,AVGPOW.powspctrm(ismember(AVGPOW.label,soimag{s}),:))
xlabel('frequency (Hz)')
ylabel('relative power change')
title('magnetometers OI')
%legend(soimag{s})
xlim([56 80])
hold on
plot(60,AVGPOW.powspctrm(ismember(AVGPOW.label,soimag{s}),AVGPOW.freq==60),'r*')
hold on
plot(67,AVGPOW.powspctrm(ismember(AVGPOW.label,soimag{s}),AVGPOW.freq==67),'r*')

subplot(212)
plot(AVGPOW.freq,AVGPOW.powspctrm(ismember(AVGPOW.label,soicmb{s}),:))
xlim([56 80])
xlabel('frequency (Hz)')
ylabel('relative power change')
title('cmb grads OI')
%legend(soicmb{s})
hold on
plot(60,AVGPOW.powspctrm(ismember(AVGPOW.label,soicmb{s}),AVGPOW.freq==60),'r*')
hold on
plot(67,AVGPOW.powspctrm(ismember(AVGPOW.label,soicmb{s}),AVGPOW.freq==67),'r*')

print(fig,fullfile(figpth,subjfolds{s},'avg_tfr_mags_cmb_raw'),'-dpng')

clear TFRall