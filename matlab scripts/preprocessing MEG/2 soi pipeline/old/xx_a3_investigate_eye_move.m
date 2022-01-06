%% VS + RFT
% PhD project 2

% inspect eye movement data

% [c] Katharina Duecker

clear all; close all; clc; beep off
% define paths
pth = 'E:\UoB\Proj2_Visual Search';
mtpth = 'Y:\Visual Search RFT\matlab scripts';
addpath(fullfile(mtpth,'edf-converter-master')); %edf2mat converter
% addpath('Y:\Visual Search RFT\matlab scripts\edfconverter')
addpath('C:\Program Files (x86)\SR Research\EyeLink\bin')
command = 'edf2asc.exe';
dtpth = fullfile(pth,'data');
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')            % fieldtrip
addpath('Y:\Visual Search RFT\matlab scripts\cbrewer')
pltpth = fullfile('Y:\Visual Search RFT\results\eyelink','figures');
elpth = fullfile('Y:\Visual Search RFT\results\eyelink','1 identify saccades');
prepropth = 'Y:\Visual Search RFT\results\meg\2 merged edf mat';

mkdir(pltpth)
mkdir(elpth)
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds
% screen settings: pixel -> degree
propixx_res = [1920 1080];                                   % propixx resolution
scr.w            = 72;                                       % screen width in cm
scr.h            = 40.5;                                     % screen height in cm
scr.d            = 142;
scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
scr.scrdegrw     = asind(scr.w/scr.ch);
scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pi

% load trigger
load(fullfile(pth,'experiment','trigdef.mat'))

% set up such that centre coordinates are centre of histogram bin
hist_bins_x = [sort(propixx_res(1)/2-scr.onedegrpix/4:-scr.onedegrpix/4:0),propixx_res(1)/2+scr.onedegrpix/4:scr.onedegrpix/4:propixx_res(1)];
hist_bins_y = [sort(propixx_res(2)/2-scr.onedegrpix/4:-scr.onedegrpix/4:0),propixx_res(2)/2+scr.onedegrpix/4:scr.onedegrpix/4:propixx_res(2)];

for s = 1%:12%:length(subjfolds)
    trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
    elpth = fullfile('Y:\Visual Search RFT\results\eyelink','1 identify saccades');

    if exist(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'))
        load(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'))
        load(fullfile(trl_merge_pth,subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'elinfo')
        
        el_fsample = elinfo.fsample;
        clear elinfo
       
    else
        try
            d = dir(fullfile(dtpth,subjfolds{s}));
            f = {d.name};
            % find edf file
            idx = cellfun(@(x) regexp(x,'.edf'),f,'UniformOutput',false);
            idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
            if ~exist(fullfile(dtpth,subjfolds{s},[f{idxx}(1:end-4),'.asc']))
                cd(fullfile(dtpth,subjfolds{s}))
                system([command ' ' f{idxx}])
                cd(mtpth)
            end
            el = Edf2Mat(fullfile(dtpth,subjfolds{s},f{idxx}));
            
             % sampling rate
        cfg = [];
        cfg.dataset = fullfile(dtpth,subjfolds{s},[f{idxx}(1:end-4),'.asc']);
        el_ft = ft_preprocessing(cfg);
        
        
        if el_ft.fsample > 900
            el_fsample = round(el_ft.fsample, -3);
        else
            el_fsample = round(el_ft.fsample, -2);
        end
        clear el_ft
        catch ME
            disp('no el file');
            continue
        end
    end
    
    % List trial type
    % extract el triggers
    idx = cellfun(@(x) regexp(x,'trig:'),el.Events.Messages.info,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    eltrials = el.Events.Messages.info(idxx);
    elsamples = el.Events.Messages.time(idxx);
    
    % extract trigger values
    extrtrig = cell2mat(cellfun(@str2num,cellfun(@(x) x(7:end),eltrials,'UniformOutput',false),'UniformOutput',false));
    
    % if practice trials are included (more than 24 1's = start of block) ->
    % only select experimental trials (last 24)
    
    begexp = find(extrtrig == 1);
    begexp = begexp(end-23:end);
    extrtrig = extrtrig(begexp(1):end);
    elsamples = elsamples(begexp(1):end);
    eltrials = eltrials(begexp(1):end);
    
    % eltrl: trigger, sample begin, sample end
    
    
    eltrl = [];
    a = 1;
    while a+2 <= length(extrtrig)
        % a: beginning of sample, a+1: onset search display, a+2: trial end (4)
        
        if extrtrig(a) >= 5 && extrtrig(a) <= trigdef{end,1} && extrtrig(a+1) == 2 && extrtrig(a+2) == 4
            eltrl = [eltrl;extrtrig(a) elsamples(a) elsamples(a+1) elsamples(a+2)];
        end
        a = a + 1;
    end
    
    eltrl = [eltrl,zeros(size(eltrl,1),3)];
    elrt = zeros(size(eltrl,1),1);
    
    fixdot = propixx_res./2; % coordinates fixation do
    elxy = zeros(size(eltrl,1),5.5.*el_fsample,2);
    for e = 1:size(eltrl,1)
        % start sample: index of time
        bs = find(el.Samples.time == eltrl(e,2));
        vs = find(el.Samples.time == eltrl(e,3));
        es = find(el.Samples.time == eltrl(e,4));
        % store sample info for trial
        eltrl(e,end-2:end) = [bs vs es];
        elxy(e,1:es-bs+1,1) = el.Samples.posX(bs:es)';
        elxy(e,1:es-bs+1,2) = el.Samples.posY(bs:es)';
        
        % check rt
        elrt(e) = (vs-bs)/el_fsample;
        
        
    end
    
%     % 1. saccades
%     sacrej = [];
%     c = 0; % counter variable
%     for sac = 1:length(el.Events.Esacc.end)
%         
%         % starts within
%         y = find((eltrl(:,3) < el.Events.Esacc.start(sac)) + (el.Events.Esacc.start(sac) < eltrl(:,4) - (0.4 * el_fsample)) == 2 );
%         
%         if ~isempty(y)
%             c = c+1;
%             
%             % start,duration of saccade
%             sac_start = el.Events.Ssacc.time(sac) - eltrl(y,2);
%             sac_dur = el.Events.Esacc.duration(sac);
%             
%             sacrej(c,:) = [y,sac_start, sac_dur];
%         end
% 
%         
%     end
%     
%        
%     % 2. find saccades >1 °
%     c = 0; % counter variable
%     sacrej_liberal = [];
%     for sac = 1:length(sacrej)
%         
%         % loop over saccades - where does eye movement deviate from centre
%         % by more than 1°?
%         
%         % cut out eye position from end of baseline to end of trial
%         bsl = 1.5*el_fsample;                                % baseline dur
%         trl_l = eltrl(sacrej(sac),end)-eltrl(sacrej(sac),end-2); % trial length
%         cur_trl = squeeze(elxy(sacrej(sac),bsl:trl_l,:));
%         % deviation x direction
%         if ~isempty(find(cur_trl(:,1) < propixx_res(1)/2 - scr.onedegrpix)) || ~isempty(find(cur_trl(:,1) > propixx_res(1)/2 + scr.onedegrpix))
%             c = c+1;
%             
%             sacrej_liberal(c,:) = sacrej(sac);
%             
%             % deviation y direction
%         elseif ~isempty(find(cur_trl(:,2) < propixx_res(2)/2 - scr.onedegrpix)) || ~isempty(find(cur_trl(:,2) > propixx_res(2)/2 + scr.onedegrpix))
%             c = c+1;
%             
%             sacrej_liberal(c,:) = sacrej(sac);
%         end
%     end
    
    % 2. find eye position further away from centre than 1°
    c = 0; % counter variable
    sacrej_liberal = [];
    for t = 1:size(elxy,1)
        
        % loop over saccades - where does eye movement deviate from centre
        % by more than 1°?
        
        % cut out eye position from end of baseline to end of trial
        bsl = 1.5*el_fsample;                                % baseline dur
        trl_l = eltrl(t,end)-eltrl(t,end-2); % trial length
        cur_trl = squeeze(elxy(t,bsl:trl_l,:));
        % deviation x direction
        if ~isempty(find(cur_trl(:,1) < propixx_res(1)/2 - scr.onedegrpix)) || ~isempty(find(cur_trl(:,1) > propixx_res(1)/2 + scr.onedegrpix))
            c = c+1;
            
            sacrej_liberal(c) = t;
            
            % deviation y direction
        elseif ~isempty(find(cur_trl(:,2) < propixx_res(2)/2 - scr.onedegrpix)) || ~isempty(find(cur_trl(:,2) > propixx_res(2)/2 + scr.onedegrpix))
            c = c+1;
            
            sacrej_liberal(c) = t;
        end
    end
    sacrej_liberal = sacrej_liberal';
    
    % blinks
    blinkrej = [];
    for sac = 1:length(el.Events.Sblink.time)
        x = find((eltrl(:,2) <= el.Events.Sblink.time(sac)) + (eltrl(:,4) >= el.Events.Sblink.time(sac)) == 2 );
        
        blinkrej = [blinkrej;x];
    end
    blinkrej = unique(blinkrej);
    
    %% How much time do people spend at fix dot?
    
    % probability density clean trials (baseline vs during search)
    clean_trl = setdiff(1:size(eltrl,1),unique([sacrej_liberal;blinkrej]));
    
    pdf_clean = zeros(size(clean_trl,2),length(hist_bins_x)-1,length(hist_bins_y)-1);
    
    % create probability density for each trial
    for cltrl = 1:length(clean_trl)
        % extract eye position
        x_coord = elxy(clean_trl(cltrl),1.5*el_fsample:end-0.4*el_fsample,1);
        y_coord = elxy(clean_trl(cltrl),1.5*el_fsample:end-0.4*el_fsample,2);
        % discard zeros (duration up to button press)
        x_coord = x_coord(x_coord ~= 0);
        y_coord = y_coord(x_coord ~= 0);
        
        % estimate probability denisty for this trial
        [pdf_clean(cltrl,:,:),Xedges,Yedges] = histcounts2(x_coord,y_coord,hist_bins_x,hist_bins_y,'Normalization','probability');
    end
    
    % repeat for saccade trials
    pdf_sac = zeros(length(unique(sacrej_liberal)),length(hist_bins_x)-1,length(hist_bins_y)-1);
    
    for sac = 1:size(sacrej_liberal)
        x_coord = elxy(sacrej_liberal(sac),1.5*el_fsample:end-0.4*el_fsample,1);
        y_coord = elxy(sacrej_liberal(sac),1.5*el_fsample:end-0.4*el_fsample,2);
        x_coord = x_coord(x_coord ~= 0);
        y_coord = y_coord(x_coord ~= 0);
        
        [pdf_sac(sac,:,:),Xedges,Yedges] = histcounts2(x_coord,y_coord,hist_bins_x,hist_bins_y,'Normalization','probability');
    end
    
   
    % colormap
    cm = cbrewer('div','RdBu',21);
    cm = cm(1:11,:);
    
    m_pdf_clean = squeeze(mean(pdf_clean,1));
    fig = figure;
    subplot(211)
    imagesc(Xedges./scr.onedegrpix, Yedges./scr.onedegrpix,m_pdf_clean)
    hold on
    % degree around fix dot
     rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix)./scr.onedegrpix, 1*2, 1*2]);
    hold on
    % search display
    rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix*5)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix*5)./scr.onedegrpix, 5*2, 5*2])
    
    xticks(0:4:scr.scrdegrw)
    yticks(0:3:scr.scrdegrh)
    xlabel('x degree')
    ylabel('y degree')
    title('no saccades')
    colormap(flipud(cm))
    caxis([0, round(max(m_pdf_clean(:)),1)])
    cb = colorbar;
    cb.Label.String = 'probability of coordinates';
    axis xy
    pbaspect([ceil(scr.scrdegrw)/ceil(scr.scrdegrh), 1,1])
    subplot(212)
    m_pdf_sac = squeeze(mean(pdf_sac,1));
    imagesc(Xedges./scr.onedegrpix, Yedges./scr.onedegrpix,m_pdf_sac)
     % degree around fix dot
     rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix)./scr.onedegrpix, 1*2, 1*2]);
    hold on
    % search display
    rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix*5)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix*5)./scr.onedegrpix, 5*2, 5*2])
    
    xlabel('x degree')
    ylabel('y degree')
    title('saccades')
    xticks(0:4:scr.scrdegrw)
    yticks(0:3:scr.scrdegrh)
    colormap(flipud(cm))
    caxis([0, round(max(m_pdf_sac(:)),1)])
    cb = colorbar;
    cb.Label.String = 'probability of coordinates';
    axis xy
    pbaspect([ceil(scr.scrdegrw)/ceil(scr.scrdegrh), 1,1])
    
    print(fig,fullfile(pltpth,[subjfolds{s},'_heatmap_Esacc']),'-dpng')
    
    
    elinfo.fsample = el_fsample;
    elinfo.eyepos = elxy;
    elinfo.trials = eltrl;
    %elinfo.sactrial = sacrej;
    elinfo.sactrial_liv = sacrej_liberal;
    elinfo.blink = blinkrej;
    elinfo.rt = elrt;
    elinfo.pdf_clean = pdf_clean;
    elinfo.pdf_sac = pdf_sac;
    
    clear el_fsample elxy eltrl elrt pdf*
    
    if ~exist(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'))
        save(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'),'el','-v7.3')
    end
    
    save(fullfile(elpth,[subjfolds{s},'_elinfo.mat']),'elinfo')
    
    clearvars  N* b* a* x* y* c* e* f* X* Y* idx* d trl* sac*
    close all
end