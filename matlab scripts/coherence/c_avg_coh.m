%% VS + RFT
% PhD project 2

% average coherence over interval 0 - median RT

clear all; close all; clc; beep off;

% define paths
pth = 'Y:\Visual Search RFT';

% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');

% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb','coh');
% coherence soi
cohsoipth = fullfile(pth,'results','meg', '5 COH hilb','soi');

% figure path
figpth = fullfile(pth,'results','meg', 'figures', 'COH hilb');

% sampling freq
fs = 1000;
% foi
freqvec = 55:75;

% list subj
d = dir(cohpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% load soi
for s = 4%11:length(subjfolds)
    load(fullfile(cohsoipth,subjfolds{s},'soi_stat.mat'))

    % load trial info
    load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'))
    
    med_rt = median(cell2mat(rspinfo.trl(:,3)));        % median RT
    
    % load coherence results
    load(fullfile(cohpth,subjfolds{s},'coh_freqw_2_T60.mat'))
    load(fullfile(cohpth,subjfolds{s},'coh_freqw_2_T67.mat'))
    
    % select 0 - median rt
    bsl_l = 1.5;        % baseline length
    ilow = bsl_l * fs;
    iup = round(bsl_l*fs + med_rt * fs);
    
    % T 60 Hz
    avg_coh60.Tgrad = mean(coh60.coh60Tgrad(:,:,ilow:iup),3); % Target grads
    avg_coh60.Tmag = mean(coh60.coh60Tmag(:,:,ilow:iup),3);   % Target mags
    avg_coh60.Dgrad = mean(coh60.coh60Dgrad(:,:,ilow:iup),3); % distractor grads
    avg_coh60.Dmag = mean(coh60.coh60Dgrad(:,:,ilow:iup),3);  % Distractor mags
    
    
    % T 67 Hz
    avg_coh67.Tgrad = mean(coh67.coh67Tgrad(:,:,ilow:iup),3); % Target grads
    avg_coh67.Tmag = mean(coh67.coh67Tmag(:,:,ilow:iup),3);   % Target mags
    avg_coh67.Dgrad = mean(coh67.coh67Dgrad(:,:,ilow:iup),3); % distractor grads
    avg_coh67.Dmag = mean(coh67.coh67Dgrad(:,:,ilow:iup),3);  % Distractor mags
    % plots
    mkdir(fullfile(figpth,subjfolds{s}))
    % gradiometers
    
    % concatenate into matrix for plotting
    all_avg_coh_grad = [avg_coh60.Tgrad(1:end-2,:); avg_coh60.Dgrad(1:end-2,:); avg_coh67.Tgrad(1:end-2,:); avg_coh67.Dgrad(1:end-2,:)];
    all_avg_coh_mag = [avg_coh60.Tmag(1:end-2,:); avg_coh60.Dmag(1:end-2,:); avg_coh67.Tmag(1:end-2,:); avg_coh67.Dmag(1:end-2,:)];
    
    fig = figure;
    set(gcf,'Position',[0 0 1500 750])
    if length(soigrad) > 12
        rows = 4;
    else
        rows = 2;
    end
    cols = ceil((size(coh60.coh60Tgrad,1)-2)/rows);
    for c = 1:size(coh60.coh60Tgrad,1)-2
        subplot(rows,cols,c)
        plot(freqvec, avg_coh60.Tgrad(c,:),'Color',[0 0.4470 0.7410])
        hold on
        plot(freqvec, avg_coh60.Dgrad(c,:), 'Color',[0.8500 0.3250 0.0980])
        hold on
        plot(freqvec, avg_coh67.Tgrad(c,:),'Color',[0.4660 0.6740 0.1880])
        hold on
        plot(freqvec, avg_coh67.Dgrad(c,:), 'Color',[0.6350 0.0780 0.1840])
        legend({'T 60 Hz','D 67 Hz', 'T 67 Hz', 'D 60 Hz'},'AutoUpdate','off')
        % mark foi
        hold on
        plot(60,avg_coh60.Tgrad(c,freqvec == 60),'*','Color',[0 0.4470 0.7410])
        hold on
        plot(67,avg_coh60.Tgrad(c,freqvec == 67),'*','Color',[0 0.4470 0.7410])
        hold on
        plot(60,avg_coh60.Dgrad(c,freqvec == 60),'*','Color',[0.8500 0.3250 0.0980])
        hold on
        plot(67,avg_coh60.Dgrad(c,freqvec == 67),'*','Color',[0.8500 0.3250 0.0980])
        hold on
        plot(60,avg_coh67.Tgrad(c,freqvec == 60),'*','Color',[0.4660 0.6740 0.1880])
        hold on
        plot(67,avg_coh67.Tgrad(c,freqvec == 67),'*','Color',[0.4660 0.6740 0.1880])
        hold on
        plot(60,avg_coh67.Dgrad(c,freqvec == 60),'*','Color',[0.6350 0.0780 0.1840])
        hold on
        plot(67,avg_coh67.Dgrad(c,freqvec == 67),'*','Color',[0.6350 0.0780 0.1840])
        ylim([0 round(max(all_avg_coh_grad(:))+0.01,2)])
        %xlim([-1.5 3])
        xlabel('frequency (Hz)')
        ylabel('coherence')
        axis xy
        title(['average COH T 60 Hz ',soigrad{c}])
    end
    print(fig,fullfile(figpth,subjfolds{s},'avg_coh_soigrad_t'),'-dpng')
    
    % magnetometers
    fig = figure;
    set(gcf,'Position',[0 0 1500 750])
    if length(soimag) > 12
        rows = 4;
    else
        rows = 2;
    end
    cols = ceil((size(coh60.coh60Tmag,1)-2)/rows);
    for c = 1:size(coh60.coh60Tmag,1)-2
        subplot(rows,cols,c)
        plot(freqvec, avg_coh60.Tmag(c,:),'Color',[0 0.4470 0.7410])
        hold on
        plot(freqvec, avg_coh60.Dmag(c,:), 'Color',[0.8500 0.3250 0.0980])
        hold on
        plot(freqvec, avg_coh67.Tmag(c,:),'Color',[0.4660 0.6740 0.1880])
        hold on
        plot(freqvec, avg_coh67.Dmag(c,:), 'Color',[0.6350 0.0780 0.1840])
        legend({'T 60 Hz','D 67 Hz', 'T 67 Hz', 'D 60 Hz'},'AutoUpdate','off')
        % mark foi
        hold on
        plot(60,avg_coh60.Tmag(c,freqvec == 60),'*','Color',[0 0.4470 0.7410])
        hold on
        plot(67,avg_coh60.Tmag(c,freqvec == 67),'*','Color',[0 0.4470 0.7410])
        hold on
        plot(60,avg_coh60.Dmag(c,freqvec == 60),'*','Color',[0.8500 0.3250 0.0980])
        hold on
        plot(67,avg_coh60.Dmag(c,freqvec == 67),'*','Color',[0.8500 0.3250 0.0980])
        hold on
        plot(60,avg_coh67.Tmag(c,freqvec == 60),'*','Color',[0.4660 0.6740 0.1880])
        hold on
        plot(67,avg_coh67.Tmag(c,freqvec == 67),'*','Color',[0.4660 0.6740 0.1880])
        hold on
        plot(60,avg_coh67.Dmag(c,freqvec == 60),'*','Color',[0.6350 0.0780 0.1840])
        hold on
        plot(67,avg_coh67.Dmag(c,freqvec == 67),'*','Color',[0.6350 0.0780 0.1840])
        ylim([0 round(max(all_avg_coh_mag(:))+0.01,2)])
        %xlim([-1.5 3])
        xlabel('frequency (Hz)')
        ylabel('coherence')
        axis xy
        title(['avg coh T 60 Hz ',soimag{c}])
    end
    print(fig,fullfile(figpth,subjfolds{s},'avg_coh_soimag_t'),'-dpng')
    
    clear coh6* psd6* csd6* avg_* rsp* alltrl*
    close all
end