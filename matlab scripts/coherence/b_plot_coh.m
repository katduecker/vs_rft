%% VS + RFT
% PhD project 2

% loads result of coherence analysis and plots coherence over time

clear all; close all; clc; beep off;
% define paths
pth = 'Y:\Visual Search RFT';

% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');

% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb','coh');
cohsoipth = fullfile(pth,'results','meg', '5 COH hilb','soi');
figpth = fullfile(pth,'results','meg', 'figures', 'COH hilb');
mkdir(figpth)

% list subj
d = dir(cohpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% load coherence results
for s = 10%1:length(subjfolds)
    load(fullfile(cohsoipth,subjfolds{s},'soi_stat.mat'))
    load(fullfile(cohpth,subjfolds{s},'coh_freqw_2_T60.mat'))
    load(fullfile(cohpth,subjfolds{s},'coh_freqw_2_T67.mat'))
    
    %% Plots
    mkdir(fullfile(figpth,subjfolds{s}))
    timevec = linspace(-1.5,2,3501);
    freqvec = 55:75;
    
    % PSD diode
    fig = figure;
    set(gcf,'Position',[0 0 1800 1400])
    subplot(221)
    imagesc(timevec,freqvec,squeeze(psd60.psdmisc60Tgrad(1,:,:)))
    %xlim([-1.5 3])
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    axis xy
    title('Target freq 60 Hz')
    subplot(222)
    imagesc(timevec,freqvec,squeeze(psd60.psdmisc60Dgrad(1,:,:)))
    %xlim([-1.5 3])
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    axis xy
    title('Distractor freq 67 Hz')
    subplot(223)
    imagesc(timevec,freqvec,squeeze(psd67.psdmisc67Tgrad(1,:,:)))
    %xlim([-1.5 3])
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    axis xy
    title('Target freq 67 Hz')
    subplot(224)
    imagesc(timevec,freqvec,squeeze(psd67.psdmisc67Dgrad(1,:,:)))
    %xlim([-1.5 3])
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    axis xy
    title('Distractor freq 60 Hz')
    print(fig,fullfile(figpth,subjfolds{s},'psd_diode'),'-dpng')
    
    % PSD grads
    fig = figure;
    set(gcf,'Position',[0 0 1800 1400])
    cols = round((size(coh60.coh60Tgrad,1)-2)/2);
    for c = 1:size(coh60.coh60Tgrad,1)-2
        subplot(2,cols,c)
        imagesc(timevec,freqvec,squeeze(psd60.psdmeg60Tgrad(c,:,:)))
        %xlim([-1.5 3])
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        axis xy
        title(['PSD T 60 Hz ',soigrad{c}])
    end
    print(fig,fullfile(figpth,subjfolds{s},'psd_soigrad_t'),'-dpng')
    
    % PSD mags
    fig = figure;
    set(gcf,'Position',[0 0 1800 1400])
    cols = round((size(coh60.coh60Tmag,1)-2)/2);
    for c = 1:size(coh60.coh60Tmag,1)-2
        subplot(2,cols,c)
        imagesc(timevec,freqvec,squeeze(psd60.psdmeg60Tmag(c,:,:)))
        %xlim([-1.5 3])
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        axis xy
        title(['PSD T 60 Hz ',soimag{c}])
    end
    print(fig,fullfile(figpth,subjfolds{s},'psd_soimag_t'),'-dpng')
    
    % COH grads
    fig = figure;
    set(gcf,'Position',[0 0 1800 1400])
    cols = round((size(coh60.coh60Tgrad,1)-2)/2);
    for c = 1:size(coh60.coh60Tgrad,1)-2
        subplot(2,cols,c)
        imagesc(timevec,freqvec,squeeze(coh60.coh60Tgrad(c,:,:)))
        %xlim([-1.5 3])
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        axis xy
        title(['COH T 60 Hz ',soigrad{c}])
    end
    print(fig,fullfile(figpth,subjfolds{s},'coh_soigrad_t'),'-dpng')
    
    % COH mags
    fig = figure;
    cols = round((size(coh60.coh60Tmag,1)-2)/2);
    for c = 1:size(coh60.coh60Tmag,1)-2
        subplot(2,cols,c)
        imagesc(timevec,freqvec,squeeze(coh60.coh60Tmag(c,:,:)))
        %xlim([-1.5 3])
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        axis xy
        title(['COH T 60 Hz ',soimag{c}])
    end
    print(fig,fullfile(figpth,subjfolds{s},'coh_soimag_t'),'-dpng')
    
    clear coh6* psd6* csd6*
    close all
end