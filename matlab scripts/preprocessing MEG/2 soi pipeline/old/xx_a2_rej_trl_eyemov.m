%% VS + RFT
% PhD project 2

% check which trials to discard based on eye movement

clear all; close all; clc; beep off
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
dtpth = fullfile(pth,'data');
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;
propixx_res = [1920 1080];              % propix

% pixel in degree
scr.w            = 72;                                       % screen width in cm
scr.h            = 40.5;                                     % screen height in cm
scr.d            = 142; 
scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pix


% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

%edf_dat_av = zeros(length(subjfolds),1);
%save('edf_dat_av.mat','edf_dat_av')
load('edf_dat_av.mat')
% 
% for s = 1:length(subjfolds)
%     d = dir(fullfile(dtpth,subjfolds{s}));
%     files = {d.name};
%     % find edf file
%     f_edf = cell2mat(cellfun(@(x) ~isempty(x),cellfun(@(x) strfind(x,'.edf'),files,'UniformOutput',false),...
%         'UniformOutput',false));
%     if find(f_edf)
%         edf_dat_av(s) = 1;
%     end
% end
% 
% % subjects to be excluded
% excl_subj = subjfolds(~edf_dat_av);   % documented in excel file
% 
% % subjects to be kept
% subjfolds = subjfolds(logical(edf_dat_av));

%% check eye movement
for s = 2%:length(subjfolds)
    load(fullfile(trl_merge_pth, subjfolds{s},'trl_overlap_meg_el_rsp_sac_search.mat'))
    
    % convert x,y pixel coordinates into degree
    
    % check eye movement vs no eye movement trials
    
    close all
    % randomly pick 5 eye movement trials
    sac_trl = elinfo.saclib(randperm(size(elinfo.saclib,1),5));
    % randomly select 5 clean trials
    clean_trl = find(~ismember(1:size(alltrl_list,1),unique([elinfo.saclib;elinfo.blink])));              % all clean trials
    clean_trl_select = clean_trl(randperm(size(clean_trl,2),5));
    timevec = linspace(-1.5, 4, elinfo.fsample*5.5);
    
    figure;
    heatmap(squeeze(mean(elinfo.move(clean_trl_select,:,:))))
    
    figure;
    subplot(221)
    for c = 1:length(clean_trl_select)
        plot(timevec,squeeze(elinfo.move(clean_trl_select(c),:,1)))
        hold on
        line(timevec,repmat(propixx_res(1)/2,1,length(timevec)) + scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        line(timevec,repmat(propixx_res(1)/2,1,length(timevec)) - scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        xlim([-1.5 3])
        ylim([0 1920])
    end
    title('clean trials x')
    subplot(222)
    for c = 1:length(clean_trl_select)
        plot(timevec,squeeze(elinfo.move(clean_trl_select(c),:,2)))
        hold on
        line(timevec,repmat(propixx_res(2)/2,1,length(timevec)) + scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        line(timevec,repmat(propixx_res(2)/2,1,length(timevec)) - scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        
        xlim([-1.5 3])
        ylim([0 1080])

    end
    title('clean trials y')
    subplot(223)
    for c = 1:length(sac_trl)
        plot(timevec,squeeze(elinfo.move(sac_trl(c),:,1)))
        hold on
        line(timevec,repmat(propixx_res(1)/2,1,length(timevec)) + scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        line(timevec,repmat(propixx_res(1)/2,1,length(timevec)) - scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        ylim([0 1920])
        xlim([-1.5 3])

    end
    title('eye movement trials x')
    subplot(224)
    for c = 1:length(sac_trl)
        plot(timevec,squeeze(elinfo.move(sac_trl(c),:,2)))
        hold on
        line(timevec,repmat(propixx_res(2)/2,1,length(timevec)) + scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        line(timevec,repmat(propixx_res(2)/2,1,length(timevec)) - scr.onedegrpix,'Color','k','LineStyle','-.')
        hold on
        ylim([0 1080])
        xlim([-1.5 3])
    end
    title('eye movement trials y')
    
end
