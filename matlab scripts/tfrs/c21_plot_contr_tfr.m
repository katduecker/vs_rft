%% VS + RFT
% PhD project 2

% plot averaged contrast TFR (Target 60 Hz vs Target 67 Hz)

% [c] Katharina Duecker

clear all; close all; clc; beep off;
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');
% save tfr
tfrpth = fullfile(pth,'results','meg', '4 TFR power');
% figure paths
figpth = fullfile(pth,'results','meg', 'figures', 'TFR power');
mkdir(tfrpth)
mkdir(figpth)
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

load('soi_tfr_subj.mat')

for s = [1:3,5:8]%length(subjfolds)
    
    %load in contrast data
    load(fullfile(tfrpth,[subjfolds{s},'_TFR_Tfreq.mat']))
    
    % compute average relative power change over each condition
    
    % apply baseline
%     cfg = [];
%     cfg.baseline = [-1 -.5];
%     cfg.baselinetype = 'relchange';
%     TFR60Tbsl = ft_freqbaseline(cfg, TFR60T);
%     TFR67Tbsl = ft_freqbaseline(cfg, TFR67T);
    
    % average over time
    cfg = [];
    cfg.latency = [.5 1];
    cfg.avgovertime = 'yes';
    TFR60Tavg = ft_selectdata(cfg, TFR60T);
    TFR67Tavg = ft_selectdata(cfg, TFR67T);
    
    % relative contrast
    TFRavgcontr_rel = TFR60Tavg;
    TFRavgcontr_rel.powspctrm = TFR60Tavg.powspctrm./TFR67Tavg.powspctrm - 1;
    
    % absolute contrast
    TFRavgcontr_abs = TFR60Tavg;
    TFRavgcontr_abs.powspctrm = TFR60Tavg.powspctrm - TFR67Tavg.powspctrm;
    
    % plot
    fig = figure;
    subplot(221)
    plot(TFRavgcontr_rel.freq,TFRavgcontr_rel.powspctrm(ismember(TFRavgcontr_rel.label,soimag{s}),:))
    xlabel('frequency (Hz)')
    ylabel('relative contrast')
    title('relative contrast magnOI')
    hold on
    plot(60,TFRavgcontr_rel.powspctrm(ismember(TFRavgcontr_rel.label,soimag{s}),TFRavgcontr_rel.freq==60),'r*')
    hold on
    plot(67,TFRavgcontr_rel.powspctrm(ismember(TFRavgcontr_rel.label,soimag{s}),TFRavgcontr_rel.freq==67),'r*')
    
    xlim([58 69])
    subplot(222)
    plot(TFRavgcontr_rel.freq,TFRavgcontr_rel.powspctrm(ismember(TFRavgcontr_rel.label,soicmb{s}),:))
    xlim([58 69])
    xlabel('frequency (Hz)')
    ylabel('relative contrast')
    title('relative contrast gradOI')
    hold on
    plot(60,TFRavgcontr_rel.powspctrm(ismember(TFRavgcontr_rel.label,soicmb{s}),TFRavgcontr_rel.freq==60),'r*')
    hold on
    plot(67,TFRavgcontr_rel.powspctrm(ismember(TFRavgcontr_rel.label,soicmb{s}),TFRavgcontr_rel.freq==67),'r*')
    
    
    subplot(223)
    plot(TFRavgcontr_abs.freq,TFRavgcontr_abs.powspctrm(ismember(TFRavgcontr_abs.label,soimag{s}),:))
    xlabel('frequency (Hz)')
    ylabel('absolute contrast')
    title('absolute contrast magnOI')
    xlim([58 69])
    hold on
    plot(60,TFRavgcontr_abs.powspctrm(ismember(TFRavgcontr_abs.label,soimag{s}),TFRavgcontr_abs.freq==60),'r*')
    hold on
    plot(67,TFRavgcontr_abs.powspctrm(ismember(TFRavgcontr_abs.label,soimag{s}),TFRavgcontr_abs.freq==67),'r*')
    
    
    subplot(224)
    plot(TFRavgcontr_abs.freq,TFRavgcontr_abs.powspctrm(ismember(TFRavgcontr_abs.label,soicmb{s}),:))
    xlim([58 69])
    xlabel('frequency (Hz)')
    ylabel('absolute contrast')
    title('absolute contrast gradOI')
    hold on
    plot(60,TFRavgcontr_abs.powspctrm(ismember(TFRavgcontr_abs.label,soicmb{s}),TFRavgcontr_abs.freq==60),'r*')
    hold on
    plot(67,TFRavgcontr_abs.powspctrm(ismember(TFRavgcontr_abs.label,soicmb{s}),TFRavgcontr_abs.freq==67),'r*')
    
    print(fig,fullfile(figpth,subjfolds{s},'contrast_avg_tfr_changed_xlim'),'-dpng')
    
    close all
end