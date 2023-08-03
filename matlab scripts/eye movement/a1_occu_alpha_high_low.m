%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a1. occular artefacts for alpha high vs low

% Input
% -s : subject index
% -time_oit: time interval for median split


% Output
% number of saccades, blinks and gaze bias for alpha high vs low

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Eye movement analysis
% a: ocular artefacts for alpha high low
% b: ocular artefacts for fast vs slow trials

function a1_occu_alpha_high_low(s,time_oi)

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
alphapth = fullfile(pth,'results','meg','6 Alpha','not align','pow align iaf');
alphasoipth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');

inpth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid');
mergepth = fullfile(pth,'results','meg','2 merged edf mat');
dtpth = fullfile(pth,'data');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions','alpha RFT');
occupth = fullfile(pth,'results','eyelink');
mkdir(occupth)
addpath(genpath(fullfile(pth,'matlab scripts')))
condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));

folds = subj;

% saccades matrix
sac_subj_condi_bslsearch_hl = zeros(length(condi),2,2);
% blink matrix
bl_subj_condi_bslsearch_hl = zeros(length(condi),2,2);

% eye bias matrix
bias_subj_condi_hl = zeros(length(condi),2);


% alpha sensors
load(fullfile(alphasoipth,folds{s},'iaf_soi_not_align.mat'))



%% load eye movement data
d = dir(fullfile(inpth,folds{s}));
d = {d.name};
files = d(strncmp(d,'dat',3));

load(fullfile(mergepth,folds{s},"trl_overlap_meg_el_rsp.mat"))

eye_move = elinfo.move(rspinfo.keeptrl_rsp,:,:);
eye_eltrl= elinfo.eltrl(rspinfo.keeptrl_rsp,:);

eye_move = eye_move(meginfo.keeptrl_all,:,:);
eye_eltrl = eye_eltrl(meginfo.keeptrl_all,:);


% sanity - are these actually the eye coordinates in pixels?? -> to compare
% % to coordinates of stimuli
% scr.w            = 72;                                       % screen width in cm
% scr.h            = 40.5;                                     % screen height in cm
% scr.d            = 142;
% scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
% scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
% scr.scrdegrw     = asind(scr.w/scr.ch);
% scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pi
% 
% % histogram
% propixx_res = [1920 1080];                                   % propixx resolution
% 
% hist_bins_x = [sort(propixx_res(1)/2-scr.onedegrpix/4:-scr.onedegrpix/4:0),propixx_res(1)/2+scr.onedegrpix/4:scr.onedegrpix/4:propixx_res(1)];
% hist_bins_y = [sort(propixx_res(2)/2-scr.onedegrpix/4:-scr.onedegrpix/4:0),propixx_res(2)/2+scr.onedegrpix/4:scr.onedegrpix/4:propixx_res(2)];
% 
% x_coord = eye_move(:,1500:2000,1);
% y_coord = eye_move(:,1500:2000,2);
% 
% for t = 1:size(x_coord,1)
%    [pdf_eye(t,:,:),Xedges,Yedges] = histcounts2(x_coord,y_coord,hist_bins_x,hist_bins_y,'Normalization','probability');
% 
% end
% 
% fig = figure;
% imagesc(Xedges, Yedges,squeeze(mean(pdf_eye,1)))

%% load response file with info on search display
d = dir(fullfile(dtpth,folds{s}));
d = {d.name};

% find matfile (there are some irregularities in how it's named, but it
% always contains subject code & .mat)
mat_file = d{(strncmp(d,folds{s}(end-3:end),4) + cell2mat(cellfun(@(x) ~isempty(x), regexp(d,'.mat'),'UniformOutput',false))) == 2};
load(fullfile(dtpth,folds{s},mat_file))


% bring trial defintions & search displays in correct order

trl_def = {};
trl_search_dsp = {};

% loop over blocks
for bl_idx = 1:size(subj.exp.trials{1},1)

    % store current trial definition
    cur_trl = cell(1,size(subj.exp.trials{1},2));
    for trl_idx = 1:size(subj.exp.trials{1},2)
        cur_trl{trl_idx} = strjoin(cellfun(@num2str, subj.exp.trials{1}{bl_idx,trl_idx}, 'UniformOutput',false),'');
    end

    trl_def = [trl_def;cur_trl'];
    % store search displays in this block
    trl_search_dsp = [trl_search_dsp;subj.srchdsp(bl_idx,:)'];
end

% decode trial info into trigger
% are triggers the same as in rspinfo?
load(fullfile(pth,'experiment','trigdef.mat'))

% find triggers related to trial definitions
trl_trig = zeros(size(trl_def));

for tg = 5:size(trigdef,1)

    idx = ismember(trl_def,trigdef{tg,2});

    trl_trig(idx) = trigdef{tg,1};
end

% same trial order as in rspinfo?
if ~isequal([rspinfo.trl{:,1}]',trl_trig)
    error('error in reshaping experiment file!')
end

% select trials to be kept
trl_trig = trl_trig(rspinfo.keeptrl_rsp);
trl_search_dsp = trl_search_dsp(rspinfo.keeptrl_rsp);

trl_trig = trl_trig(meginfo.keeptrl_all);
trl_search_dsp = trl_search_dsp(meginfo.keeptrl_all);

trl_trig = num2cell(trl_trig);

trl_search_dsp = [trl_trig,trl_search_dsp];

clear trl_def trl_trig
%% get saccade and blink timing
load(fullfile(dtpth,folds{s},'el_struct.mat'))

sac_time = el.Events.Ssacc.time;
blink_time = el.Events.Sblink.time;

clear el

%% select condition
load(fullfile(pth, 'experiment','trigdef.mat'))


d = dir(fullfile(inpth,folds{s}));
files = {d.name};
files(1:2) = [];

for c_idx = 1:length(condi)

    % find relevant data files
    condi_files = zeros(length(files),1);
    for c = 1:length(condi{c_idx})

        condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c_idx}{c}),'UniformOutput',false))';

    end

    condi_files = condi_files == length(condi{c_idx});

    % load files 6067
    c_files = files(condi_files);

    load(fullfile(inpth,folds{s},c_files{1}));
    trl_idx = find(trlcur);

    data = data_trig;
    perf_TFR_coh = perf_cur;

    eye_move_condi = eye_move(trlcur,:,:);
    eye_trl_condi = eye_eltrl(trlcur,:,:);
    srch_dsp_condi = trl_search_dsp(trlcur,:);

    % load & append performance and eye movement
    for f = 2:length(c_files)
        clear data_trig trlcur perf_cur
        load(fullfile(inpth,folds{s},c_files{f}));

        data = ft_appenddata([],data,data_trig);
        perf_TFR_coh = [perf_TFR_coh;perf_cur];

        eye_move_condi = [eye_move_condi;eye_move(trlcur,:,:)];
        eye_trl_condi = [eye_trl_condi;eye_eltrl(trlcur,:)]; % eltrl contains samples and time points for baseline, search onset, search
        srch_dsp_condi = [srch_dsp_condi;trl_search_dsp(trlcur,:)];

    end

    % compare if data triggers and eye data triggers are the same
    if ~isequal([perf_TFR_coh{:,1}]',eye_trl_condi(:,1),[srch_dsp_condi{:,1}]')
        error('MEG and eyelink trials are not the same!')
    end


    % get alpha power for these data
    winl=0.5;
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.channel = 'MEG';
    cfg.taper = 'hanning';
    cfg.foi = 4:1/winl:30;
    cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
    cfg.toi = -1.75:0.05:1;
    cfg.keeptrials = 'yes';
    TFR_alpha = ft_freqanalysis(cfg,data);

    cfg = [];
    cfg.method = 'sum';
    TFR_alpha = ft_combineplanar(cfg,TFR_alpha);

    % average over baseline and extract IAF power
    cfg = [];
    cfg.frequency = [iaf_grad iaf_grad];
    cfg.avgoverfreq = 'yes';

    cfg.latency = time_oi;
    cfg.avgovertime = 'yes';

    cfg.channel = soi_grad_cmb;

    cfg.avgoverchan = 'yes';
    IAFpow = ft_selectdata(cfg,TFR_alpha);

    m_iaf = median(IAFpow.powspctrm);

    iaf_pow = IAFpow.powspctrm;


    trl_high = find(iaf_pow>m_iaf);
    trl_low = find(iaf_pow<m_iaf);


    %% select eye movement data for current condition

    % eye movement data and trial info for high vs low (first 500 ms)
    eye_move_high = eye_move_condi(trl_high,1500:2000-1,:);
    eye_move_low = eye_move_condi(trl_low,1500:2000-1,:);

    eye_trl_high = eye_trl_condi(trl_high,:);
    eye_trl_low = eye_trl_condi(trl_low,:);

    %% select search display for current condition

    srch_dsp_high = srch_dsp_condi(trl_high,:);
    srch_dsp_low = srch_dsp_condi(trl_low,:);

    %% find saccades that happened during trial
    s_high = zeros(size(eye_trl_high,1),2);
    s_low = zeros(size(eye_trl_low,1),2);

    bl_high = s_high;
    bl_low = s_high;

    for t = 1:size(eye_trl_high,1)

        % saccades during baseline
        s_high(t,1) = sum((eye_trl_high(t,2) < sac_time) + (eye_trl_high(t,3) > sac_time) == 2);
        s_low(t,1) = sum((eye_trl_low(t,2) < sac_time) + (eye_trl_low(t,3) > sac_time) == 2);

        % saccades during search (first 500 ms)
        s_high(t,2) = sum((eye_trl_high(t,3) < sac_time) + (eye_trl_high(t,3)+500 > sac_time) == 2);
        s_low(t,2) = sum((eye_trl_low(t,3) < sac_time) + (eye_trl_low(t,3)+500 > sac_time) == 2);


        % blinks during baseline
        bl_high(t,1) = sum((eye_trl_high(t,2) < blink_time) + (eye_trl_high(t,3) > blink_time) == 2);
        bl_low(t,1) = sum((eye_trl_low(t,2) < blink_time) + (eye_trl_low(t,3) > blink_time) == 2);

        % saccades during search
        bl_high(t,2) = sum((eye_trl_high(t,3) < blink_time) + (eye_trl_high(t,3)+500 > blink_time) == 2);
        bl_low(t,2) = sum((eye_trl_low(t,3) < blink_time) + (eye_trl_low(t,3)+500 > blink_time) == 2);


    end

    % high baseline
    sac_subj_condi_bslsearch_hl(c_idx,1,1) = mean(s_high(:,1));
    bl_subj_condi_bslsearch_hl(c_idx,1,1) = mean(bl_high(:,1));
    % high search
    sac_subj_condi_bslsearch_hl(c_idx,2,1) = mean(s_high(:,2));
    bl_subj_condi_bslsearch_hl(c_idx,2,1) = mean(bl_high(:,2));

    % low baseline
    sac_subj_condi_bslsearch_hl(c_idx,1,2) = mean(s_low(:,1));
    bl_subj_condi_bslsearch_hl(c_idx,1,2) = mean(bl_low(:,1));
    % low search
    sac_subj_condi_bslsearch_hl(c_idx,2,2) = mean(s_low(:,2));
    bl_subj_condi_bslsearch_hl(c_idx,2,2) = mean(bl_low(:,2));

    clear s_high s_low bl_high bl_low


    %% Gaze bias towards Target colour


    % bin eye movement data
    % average location in 50 ms bins
    numbin = size(eye_move_low,2)/50;

    eye_move_high_bin = zeros(length(trl_high),numbin,size(eye_move_high,3));
    eye_move_low_bin = eye_move_high_bin;
    sb = 1;     % start of bin
    for b = 1:numbin
        eye_move_high_bin(:,b,:) = mean(eye_move_high(:,sb:sb+50-1,:),2);
        eye_move_low_bin(:,b,:) = mean(eye_move_low(:,sb:sb+50-1,:),2);
        sb = sb+50;
    end

    for t = 1:size(eye_trl_high,1)



        %% Alpha high
        search_coord = srch_dsp_high{t,2}.cxy;

        % x-distance to each stimulus (x) at each time bin (y)
        x_dist = squeeze(eye_move_high_bin(t,:,1)) - search_coord(1,:)';
        % y distance
        y_dist = squeeze(eye_move_high_bin(t,:,2)) - search_coord(2,:)';

        % euclidean distance
        dist = sqrt(x_dist.^2 + y_dist.^2);

        [min_dist, p_dist] = min(dist,[],1);


        % how many stimuli where looked at and how often?
        [count_time, which_stim] = groupcounts(p_dist');

        % stimulus that eye was closest to
        [~, p] = max(count_time);
        stim_idx = which_stim(p);

        % percentage of time spent close to a stimulus in T color
        % multiply logical array (Target yes no) with number of bins
        % for each target


        stim_id_high(t) = sum(strcmp(srch_dsp_high{t,2}.stimidentity(which_stim),'lt').*(count_time'./sum(count_time)));


        %% Alpha low
        search_coord = srch_dsp_low{t,2}.cxy;

        % x-distance to each stimulus (x) at each time bin (y)
        x_dist = squeeze(eye_move_low_bin(t,:,1)) - search_coord(1,:)';
        % y distance
        y_dist = squeeze(eye_move_low_bin(t,:,2)) - search_coord(2,:)';

        % euclidean distance
        dist = sqrt(x_dist.^2 + y_dist.^2);

        [min_dist, p_dist] = min(dist,[],1);


        % how many stimuli where looked at and how often?
        [count_time, which_stim] = groupcounts(p_dist');

        % stimulus that eye was closest to
        [~, p] = max(count_time);
        stim_idx = which_stim(p);

        % percentage of time spent close to a stimulus in T color
        % multiply logical array (Target yes no) with number of bins
        % for each target


        stim_id_low(t) = sum(strcmp(srch_dsp_low{t,2}.stimidentity(which_stim),'lt').*(count_time'./sum(count_time)));

    end


    % high alpha
    bias_subj_condi_hl(c_idx,1) = mean(stim_id_high);

    % low alpha
    bias_subj_condi_hl(c_idx,2) = mean(stim_id_low);
end


time_oi_str = arrayfun(@num2str, time_oi.*1000,'UniformOutput',false);

mkdir(fullfile(occupth,folds{s}))
save(fullfile(occupth,folds{s},['alpha_high_low_',strjoin(time_oi_str,'_'),'.mat']),'sac_subj_condi_bslsearch_hl','bl_subj_condi_bslsearch_hl','bias_subj_condi_hl')
