%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a2. plot coherence for each subject (-> Supplementary figure 1);
% grandaverage (Fig 3a)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023


%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)



%% settings
clear all; close all; clc; beep off;

pth = 'Z:\Visual Search RFT';

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'SNR');

cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig', '0 SNR');
mkdir(cohfigpth)
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
datpth = fullfile(pth,'data');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath(fullfile('Z:','fieldtrip'))
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);

ft_defaults;

toi = [-2.5 2];                             % start and end of trial in sec
avgtoi = 0.5;
fw = 2;                                     % bandwidth bp filter
fs = 1000;
foi = 50:75;
% list subjects
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));


%% load in example subject to get ft structure
d = dir(fullfile(datpth,subj{1},'meg'));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(datpth,subj{1},'meg',f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = {'MEG','MISC004','MISC005'};
% load in data for this part
data = ft_preprocessing(cfg);

% fourier transform
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
%cfg.pad        = 'nextpow2';
cfg.taper      = 'hanning';
cfg.toi        = -1.5:2;
cfg.foi        = foi;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*.1;
cfg.tapsmofrq  = 3;
cfg.keeptrials = 'no';
cfg.channel    = {'MEG', 'MISC004','MISC005'};
freq       = ft_freqanalysis(cfg, data);

cfg = [];
cfg.channel = {'MEGGRAD', 'MISC004','MISC005'};
freqgrad = ft_selectdata(cfg,freq);
cfg.channel = {'MEGMAG', 'MISC004','MISC005'};
freqmag        = ft_selectdata(cfg, data);
freqgrad.time  = toi(1):1/fs:toi(2);
freqgrad.freq  = foi;
freqmag.time  = toi(1):1/fs:toi(2);
freqmag.freq  = foi;
freqmag.dimord = freqgrad.dimord;
clear data

%% Grads

% Topos
cfg = [];
cfg.marker = 'off';
cfg.zlim = 'zeromax';
cfg.comment = 'no';
cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
cfg.colorbar = 'yes';
cfg.ylim = [60 60];

coh_all = cell(1,length(subj));
for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},'SNR_MEGGRAD.mat'),'coh60')                    % load coherence
    %     load(fullfile(mergepth,subj{s},'trl_overlap_meg_el_rsp.mat'))               % load trial info
    %
    %     % minimum reaction time for plotting
    %     rt_subj = [rspinfo.trl{:,3}];
    %     rt_subj = rt_subj(rspinfo.keeptrl_rsp);
    %     rt_subj = rt_subj(meginfo.keeptrl_all);
    %     minrt(s) = min(rt_subj);
    coh_all{s} = coh60;
    clear coh60
end

coh_for_plot = freqgrad;


% This is Supplementary Fig. 1
fig = figure('Position',[0 0 920 1080]);
for s = 1:length(subj)
        load(fullfile(soipth,subj{s},'soi_stat_not_align.mat'))

    % paste into freq structure and plot
    coh_for_plot.powspctrm = coh_all{s}.*10^2;

    cfg.layout = 'neuromag306planar_helmet.mat';
    cfg.xlim = [0.0 0.5];
    cfg.highlight = 'on';
    cfg.highlightsymbol = 'o';
    cfg.highlightcolor = [5,113,176]./255;
    cfg.highlightsize = 4;
    cfg.highlightchannel = find(ismember(coh_for_plot.label,soigrad));
    subplot(8,4,s)
    ft_topoplotTFR(cfg,coh_for_plot)
    cb = colorbar;
    cb.Ticks = cb.Limits;
    cb.TickLabels = strsplit(sprintf('%0.1f ',cb.Ticks));

end

print(fig,fullfile(cohfigpth,'topo_coh60_allsubj_grad_fwdth5'),'-dpng','-r300')
print(fig,fullfile(cohfigpth,'topo_coh60_allsubj_grad_fwdth5'),'-dsvg','-r600','-painters')

close all

% this is the same for RIFT at 67, but it doesn't look good
fig = figure('Position',[0 0 920 1080]);
for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},'SNR_MEGGRAD.mat'),'coh67')                    % load coherence
    load(fullfile(mergepth,subj{s},'trl_overlap_meg_el_rsp.mat'))       % load trial info
   
    % minimum reaction time for plotting
    rt_subj = [rspinfo.trl{:,3}];
    rt_subj = rt_subj(rspinfo.keeptrl_rsp);
    rt_subj = rt_subj(meginfo.keeptrl_all);
    
    % paste into freq structure and plot
    coh_for_plot = freqgrad;
    coh_for_plot.powspctrm = coh67;

    cfg.layout = 'neuromag306planar_helmet.mat';
    cfg.xlim = [0 minrt(s)];

    subplot(9,4,s)
    ft_topoplotTFR(cfg,coh_for_plot)
    
    clear coh67 rt_subj coh_for_plot
end

print(fig,fullfile(cohfigpth,'topo_coh67_allsubj_grad_fwdth5'),'-dpng','-r300')
print(fig,fullfile(cohfigpth,'topo_coh67_allsubj_grad_fwdth5'),'-dsvg','-r600')

close all


% coherence spectra (not in manuscript)

fig = figure('Position',[0 0 920 1080]);
for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},'SNR_MEGGRAD.mat'))

    load(fullfile(mergepth,subj{s},'trl_overlap_meg_el_rsp.mat'))       % load trial info
    
    % load SOI
    load(fullfile(soipth,subj{s},'soi_stat.mat'))
    soimag = soi_stat(logical(cell2mat(cellfun(@(x) strcmp(x(end),'1'),soi_stat,'UniformOutput',false))));
    soigrad = soi_stat(~ismember(soi_stat,soimag));
    
    [~,startsamp] = min(abs(freqgrad.time - 0.15));
    [~,endsamp] = min(abs(freqgrad.time - minrt(s)));

    
    % paste into freq structure and plot
    coh_toi = mean(coh60(:,:,startsamp:endsamp),3,'omitnan');
    coh_soi = mean(coh_toi(ismember(freqgrad.label,soigrad),:),1,'omitnan');

    subplot(9,4,s)
    plot(freqgrad.freq,coh_soi)
    hold on
    
    % paste into freq structure and plot
    coh_toi = mean(coh67(:,:,startsamp:endsamp),3,'omitnan');
    coh_soi = mean(coh_toi(ismember(freqgrad.label,soigrad),:),1,'omitnan');

    subplot(9,4,s)
    plot(freqgrad.freq,coh_soi)
    clear coh60 coh67 rt_subj coh_for_plot
end

print(fig,fullfile(cohfigpth,'coh_allsubj_grad_spect_fwdth5'),'-dpng','-r300')
print(fig,fullfile(cohfigpth,'coh_allsubj_grad_spect_fwdth5'),'-dsvg','-r600')

close all

%% Grandavergae (Fig 3a)

all_subj_60 = cell(1,length(subj));
all_subj_67 = cell(1,length(subj));

cfg = [];
for s = 1:length(subj)

    load(fullfile(cohpth,subj{s},'SNR_MEGGRAD.mat'))

    freqgrad.powspctrm = coh60;
    all_subj_60{s} = freqgrad;

end

GA = ft_freqgrandaverage([],all_subj_60{:});

fig = figure;
cfg = [];
cfg.marker = 'off';
% cfg.highlight = 'labels';
% cfg.highlightsymbol = 'o';
% cfg.highlightsize = 4;
% cfg.highlightcolor = [0 0.5 0];
cfg.zlim = 'zeromax';
cfg.comment = 'no';
cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
cfg.colorbar = 'yes';
cfg.ylim = [60 60];
cfg.xlim = [0 mean(minrt)];
cfg.zlim = 'maxmin';
cfg.layout = 'neuromag306planar_helmet.mat';
ft_topoplotTFR(cfg,GA)   
cb = colorbar;
cb.Ticks = cb.Limits;
cb.TickLabels = strsplit(sprintf('%0.3f ',cb.Ticks));

print(fig,fullfile(cohfigpth,'topo_GA_grad_60_fwdth5'),'-dpng','-r300')
print(fig,fullfile(cohfigpth,'topo_GA_grad_60_fwdth5'),'-dsvg','-r0')

