%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d3. Cluster-based permutation t test, comparing coherence in alpha high
% vs low (including plots) -> Fig. 4

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023
 

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)

clear all; close all; clc


split_ta_tp = 0;            % split for target present/absent?
toi_alpha_split = [-1 0];   % alpha toi
bslcor = 1;                 % baseline corrected coherence?

if bslcor
    bslcor_suf = 'bslcor';
else
    bslcor_suf = '';
end

pth = 'Z:\Visual Search RFT';
addpath(fullfile(pth,'Violinplot-Matlab-master'))
addpath(fullfile('Z:','fieldtrip'))
ft_defaults;

set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Arial')
col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','alpha high low');

cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');
col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;
toi_alpha_split = arrayfun(@num2str,toi_alpha_split.*1000,'UniformOutput',false);

if split_ta_tp

    condi = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},{'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};
    split_suf = '_ta_tp';
    varnames = {'ungui 16', 'gui 16', 'ungui 32', 'gui 32'};

    
else
    condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
    split_suf = '';

    varnames = {'ungui 16', 'gui 16', 'ungui 32', 'gui 32'};
end

load(fullfile(cohpth,['coh_alpha_',bslcor_suf,'_',strjoin(toi_alpha_split,'_'),split_suf,'.mat']))


% time vector
min_rt = length(coh_condiT{1,2})/1000-0.5;
timevec = linspace(-0.5,min_rt,length(coh_condiT{1,2}));

% gradnaverage
for c = 1:length(condi)
    coh_condiT{c,4} = mean(coh_condiT{c,2},1);
    coh_condiT{c,5} = mean(coh_condiT{c,3},1);

    coh_condiD{c,4} = mean(coh_condiD{c,2},1);
    coh_condiD{c,5} = mean(coh_condiD{c,3},1);
end


%% Cluster-based permutation t test

% load in example subject to get ft structure
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
fs = 1000;
% list subjects
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

d = dir(fullfile(maxfpth,subj{1}));
fast = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),fast,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
fast = fast(idxx);
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(maxfpth,subj{1},fast{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = 'MEGGRAD';
% load in data for this part
data = ft_preprocessing(cfg);

% timelock analysis: get ERF strutcure

ERF = ft_timelockanalysis([],data);
ERF.time = timevec;
ERF.avg = repmat(timevec,204,1);
ERF.var = ERF.var(1:204,1:length(timevec));
ERF.dof = ERF.dof(1:204,1:length(timevec));
clear d fast idx idxx trlstruct

% for set size 32
Th = cell(1,length(subj));
Dh = cell(1,length(subj));
Tl = cell(1,length(subj));
Dl = cell(1,length(subj));


Thvslow = cell(1,length(subj));
Dhvslow = cell(1,length(subj));

if split_ta_tp
    fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

        timevec_stats = timevec(500:end);

        tcl = tiledlayout(2,2);
    for c = 1:length(condi)/2

    for s = 1:length(subj)


        condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi{c},'_'));

        % high
        T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));

        ERF.avg(1,:) = T;
        Th{s} = ERF;

        ERF.avg(1,:) = D;
        Dh{s} = ERF;


        % low
        T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));

        ERF.avg(1,:) = T;
        Tl{s} = ERF;

        ERF.avg(1,:) = D;
        Dl{s} = ERF;

        % high vs low
        cT = abs(squeeze(coh_condiT{condi_idx_ti,2}(s,:))-squeeze(coh_condiT{condi_idx_ti,3}(s,:)));
        ERF.avg(1,:) = cT;
        Thvslow{s} = ERF;

        cD = abs(squeeze(coh_condiD{condi_idx_ti,2}(s,:))-squeeze(coh_condiD{condi_idx_ti,3}(s,:)));
        ERF.avg(1,:) = cD;
        Dhvslow{s} = ERF;


    end

    % T-test Target vs Distractor

    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.minnbchan        = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 5000;
    cfg.latency = [0.0 min_rt];
    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare Target and Distractor fast vs slow

    cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail =-1;

    [statT{c}] = ft_timelockstatistics(cfg, Th{:}, Tl{:});

    [statD{c}] = ft_timelockstatistics(cfg, Dh{:}, Dl{:});

    nexttile
    plot(timevec,coh_condiT{condi_idx_ti,4},'Color',col_palette(1,:),'LineWidth',3)
    hold on
    plot(timevec,coh_condiT{condi_idx_ti,5},'Color',col_palette(2,:),'LineWidth',3)
    plot(timevec,coh_condiD{condi_idx_ti,4},'Color',col_palette(3,:),'LineWidth',3)
    plot(timevec,coh_condiD{condi_idx_ti,5},'Color',col_palette(4,:),'LineWidth',3)


    plot(timevec_stats(statT{c}.mask),statT{c}.mask(statT{c}.mask).*0.0275,'Color',col_palette(1,:),'LineWidth',5)
    plot(timevec_stats(statD{c}.mask),statD{c}.mask(statD{c}.mask).*0.03,'Color',col_palette(3,:),'LineWidth',5)

    title(varnames{c})
    end

print(fig,fullfile(cohfigpth,['rift_alpha_cluster_absent',strjoin(toi_alpha_split,'_'),bslcor_suf]),'-dpng')

    title(tcl,'target absent')

    fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

    tcl = tiledlayout(2,2);
    for c = 1:length(condi)/2

    for s = 1:length(subj)


        condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi{c+4},'_'));

        % high
        T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));

        ERF.avg(1,:) = T;
        Th{s} = ERF;

        ERF.avg(1,:) = D;
        Dh{s} = ERF;


        % low
        T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));

        ERF.avg(1,:) = T;
        Tl{s} = ERF;

        ERF.avg(1,:) = D;
        Dl{s} = ERF;

        % high vs low
        cT = abs(squeeze(coh_condiT{condi_idx_ti,2}(s,:))-squeeze(coh_condiT{condi_idx_ti,3}(s,:)));
        ERF.avg(1,:) = cT;
        Thvslow{s} = ERF;

        cD = abs(squeeze(coh_condiD{condi_idx_ti,2}(s,:))-squeeze(coh_condiD{condi_idx_ti,3}(s,:)));
        ERF.avg(1,:) = cD;
        Dhvslow{s} = ERF;


    end

    % T-test Target vs Distractor

    cfg         = [];

    cfg.method           = 'montecarlo';
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.minnbchan        = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 5000;
    cfg.latency = [0.0 min_rt];
    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare Target and Distractor fast vs slow

    cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail =-1;

    [statT{c}] = ft_timelockstatistics(cfg, Th{:}, Tl{:});

    [statD{c}] = ft_timelockstatistics(cfg, Dh{:}, Dl{:});

    nexttile
    plot(timevec,coh_condiT{condi_idx_ti,4},'Color',col_palette(1,:),'LineWidth',3)
    hold on
    plot(timevec,coh_condiT{condi_idx_ti,5},'Color',col_palette(2,:),'LineWidth',3)
    plot(timevec,coh_condiD{condi_idx_ti,4},'Color',col_palette(3,:),'LineWidth',3)
    plot(timevec,coh_condiD{condi_idx_ti,5},'Color',col_palette(4,:),'LineWidth',3)


    plot(timevec_stats(statT{c}.mask),statT{c}.mask(statT{c}.mask).*0.0275,'Color',col_palette(1,:),'LineWidth',5)
    plot(timevec_stats(statD{c}.mask),statD{c}.mask(statD{c}.mask).*0.03,'Color',col_palette(3,:),'LineWidth',5)

    title(varnames{c})

    end

    title(tcl,'target present')
    print(fig,fullfile(cohfigpth,['rift_alpha_cluster_present',strjoin(toi_alpha_split,'_'),bslcor_suf]),'-dpng')

    

else

fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

for c = 1:length(condi)

    for s = 1:length(subj)


        condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi{c},'_'));

        % fast
        T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));

        ERF.avg(1,:) = T;
        Th{s} = ERF;

        ERF.avg(1,:) = D;
        Dh{s} = ERF;


        % slow
        T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));

        ERF.avg(1,:) = T;
        Tl{s} = ERF;

        ERF.avg(1,:) = D;
        Dl{s} = ERF;

        % high vs low
        cT = abs(squeeze(coh_condiT{condi_idx_ti,2}(s,:))-squeeze(coh_condiT{condi_idx_ti,3}(s,:)));
        ERF.avg(1,:) = cT;
        Thvslow{s} = ERF;

        cD = abs(squeeze(coh_condiD{condi_idx_ti,2}(s,:))-squeeze(coh_condiD{condi_idx_ti,3}(s,:)));
        ERF.avg(1,:) = cD;
        Dhvslow{s} = ERF;


    end

    % T-test Target vs Distractor

    cfg         = [];

    cfg.method           = 'montecarlo';
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.minnbchan        = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 5000;
    cfg.latency = [0.0 min_rt];
    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare Target and Distractor fast vs slow

    cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail =-1;

    [statT{c}] = ft_timelockstatistics(cfg, Th{:}, Tl{:});

    [statD{c}] = ft_timelockstatistics(cfg, Dh{:}, Dl{:});


    cfg.tail         = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail  = -1;

    [statcontrhighlow{c}] = ft_timelockstatistics(cfg, Thvslow{:}, Dhvslow{:});


    % cohens d
    cfg         = [];

    cfg.method           = 'analytic';                                            % cluster-based
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'cohensd';                                           % within subject
    %cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution
    cfg.minnbchan        = 0;
    % cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;
    cfg.latency = [-0.15 min_rt];

    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare coherence for high and low alpha

    % Target
    cfg.clustertail = -1;
    cfg.tail        = -1;                                                       % two-sided test for both clustering and cluster-level statistic

    [statTvsDh{c}] = ft_timelockstatistics(cfg, Th{:}, Tl{:});


    % cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    % cfg.clustertail = -1;

    [statTvsDl{c}] = ft_timelockstatistics(cfg, Dh{:}, Dl{:});


    [statThlvsDhl_coh{c}] = ft_timelockstatistics(cfg, Thvslow{:}, Dhvslow{:});

       %% new plot


    timevec_stats = timevec(500:end);

    % set size 16
    subplot(2,2,c);
    plot(timevec,coh_condiT{condi_idx_ti,4},'Color',col_palette(1,:),'LineWidth',3)
    hold on
    plot(timevec,coh_condiT{condi_idx_ti,5},'Color',col_palette(2,:),'LineWidth',3)
    plot(timevec,coh_condiD{condi_idx_ti,4},'Color',col_palette(3,:),'LineWidth',3)
    plot(timevec,coh_condiD{condi_idx_ti,5},'Color',col_palette(4,:),'LineWidth',3)


    plot(timevec_stats(statT{c}.mask),statT{c}.mask(statT{c}.mask).*0.0275,'Color',col_palette(1,:),'LineWidth',5)
    plot(timevec_stats(statD{c}.mask),statD{c}.mask(statD{c}.mask).*0.03,'Color',col_palette(3,:),'LineWidth',5)

    %plot(timevec_stats(statF{c}.mask),statF{c}.mask(statF{c}.mask).*0.07,'Color',col_palette(5,:),'LineWidth',2)

    ylabel('coherence')
    xlabel('time (s)')
    xlim([-0.2 0.5])
    xticks([0:0.5:0.5])
    ylim([-0.01 0.03])
    yticks([0.0:0.03:0.03])
    title(varnames{c})
    box off

end

print(fig,fullfile(cohfigpth,['rift_alpha_cluster_',strjoin(toi_alpha_split,'_'),bslcor_suf,split_suf]),'-dpng')
print(fig,fullfile(cohfigpth,['rift_alpha_cluster_',strjoin(toi_alpha_split,'_'),bslcor_suf,split_suf]),'-dsvg')

fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

for c = 1:length(condi)
 % set size 16
    subplot(2,2,c);
    plot(statTvsDh{c}.time,statTvsDh{c}.cohensd,'Color',col_palette(1,:),'LineWidth',3,'LineStyle','-')
    hold on
    plot(statTvsDl{c}.time,statTvsDl{c}.cohensd,'Color',col_palette(4,:),'LineWidth',3,'LineStyle','-')
    

    plot(timevec_stats(statcontrhighlow{c}.mask),statcontrhighlow{c}.mask(statcontrhighlow{c}.mask).*0.7,'Color',col_palette(1,:),'LineWidth',5)

    ylabel('delta RIFT Cohen''s d')
    xlabel('time (s)')
    xlim([-0.15 0.5])
    xticks([0:0.5:0.5])
    box off
    title(varnames{c})
    ylim([-0.7 0.7])
    yticks([-0.7:0.7:0.7])
    %legend('Target alpha high vs low','Distractor alpha high vs low')
end

print(fig,fullfile(cohfigpth,['rift_alpha_cohensd_',strjoin(toi_alpha_split,'_'),bslcor_suf,split_suf]),'-dpng')
print(fig,fullfile(cohfigpth,['rift_alpha_cohensd_',strjoin(toi_alpha_split,'_'),bslcor_suf,split_suf]),'-dsvg')
% 


%% unguided average

condi = {{'ni','16t'},{'ni','32t'}};
varnames = {'ung16','ung32'};
fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

for c = 1:length(condi)

    for s = 1:length(subj)


        condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi{c},'_'));

        % fast
        T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));

        ERF.avg(1,:) = (T+D)./2;
        f_ung{s} = ERF;

        % slow
        T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));

        ERF.avg(1,:) = (T+D)./2;
        s_ung{s} = ERF;


    end

    % T-test Target vs Distractor

    cfg         = [];

    cfg.method           = 'montecarlo';
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.minnbchan        = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 5000;
    cfg.latency = [0.0 min_rt];
    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare Target and Distractor fast vs slow

    cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail =-1;

    [statf{c}] = ft_timelockstatistics(cfg, f_ung{:}, s_ung{:});




       %% new plot


    timevec_stats = timevec(650:end);

    % set size 16
    subplot(2,2,c);
    plot(timevec,(coh_condiT{condi_idx_ti,4}+coh_condiD{condi_idx_ti,4})./2,'Color',col_palette(5,:),'LineWidth',3)
    hold on
    plot(timevec,(coh_condiT{condi_idx_ti,5}+coh_condiD{condi_idx_ti,5})./2,'Color',col_palette(6,:),'LineWidth',3)

    plot(timevec_stats(statf{c}.mask),statf{c}.mask(statf{c}.mask).*0.0275,'Color',col_palette(5,:),'LineWidth',5)


    ylabel('coherence')
    xlabel('time (s)')
    xlim([-0.2 0.5])
    xticks([0:0.5:0.5])
    ylim([-0.01 0.03])
    yticks([0.0:0.03:0.03])
    title(varnames{c})
    box off

end

print(fig,fullfile(cohfigpth,['rift_alpha_cluster_ung',strjoin(toi_alpha_split,'_'),bslcor_suf,split_suf]),'-dpng')
print(fig,fullfile(cohfigpth,['rift_alpha_cluster_ung',strjoin(toi_alpha_split,'_'),bslcor_suf,split_suf]),'-dsvg')

end
% fast_TD = cell(1,length(subj));
% slow_TD = cell(1,length(subj));
% fig = figure('Position',[0 0  1980/1.5 1080/1.5]);
% 
% for c = 1:length(condi_ti)
% 
%     for s = 1:length(subj)
% 
% 
%         condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi_ti{c},'_'));
% 
%         % fast
%         T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
%         D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));
% 
%         ERF.avg(1,:) = T-D;
%         fast_TD{s} = ERF;
% 
%         coh_diff_fast(s,:) = T-D;
% 
%         % slow
%         T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
%         D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));
% 
%         ERF.avg(1,:) = T-D;
%         slow_TD{s} = ERF;
%         coh_diff_slow(s,:) = T-D;
% 
%     end
% 
%     % T-test Target vs Distractor
% 
%     cfg         = [];
% 
%     cfg.method           = 'montecarlo';                                            % cluster-based
%     cfg.channel          = ERF.label{1};
%     cfg.statistic        = 'depsamplesT';                                           % within subject
%     cfg.correctm         = 'cluster';
%     cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
%     cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution
%     cfg.minnbchan        = 0;
%     cfg.alpha            = 0.05;
%     cfg.numrandomization = 1000;
%     cfg.latency = [0 min_rt];
%     Nsubj  = length(subj);
%     design = zeros(2, Nsubj*2);
%     design(1,:) = [1:Nsubj 1:Nsubj];
%     design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
% 
%     cfg.design = design;
%     cfg.uvar   = 1;
%     cfg.ivar   = 2;
% 
%     % compare Target and Distractor fast vs slow
% 
%     cfg.tail             = 0;                                                       % two-sided test for both clustering and cluter-level statistic
%     cfg.clustertail = 0;
% 
%     [stat_fast_slow{c}] = ft_timelockstatistics(cfg, fast_TD{:}, slow_TD{:});
% 
%     subplot(2,2,c)
%     plot(timevec,mean(coh_diff_fast,1))
%     hold on 
%     plot(timevec,mean(coh_diff_slow,1))
% end

% close all
% 
% scatter(mean(coh_diff_fast(:,500:end),2),mean(coh_diff_slow(:,500:end),2))
% hold on
% plot(-0.03:0.001:0.03,-0.03:0.001:0.03)

% 
% % set size 32
% ax2 = subplot(132);
% plot(t_vec,GA_ti_32_fast(1,:),'Color',col_palette(1,:),'LineWidth',3)
% hold on
% plot(t_vec,GA_ti_32_slow(1,:),'Color',col_palette(2,:),'LineWidth',3)
% plot(t_vec,GA_ti_32_fast(2,:),'Color',col_palette(3,:),'LineWidth',3)
% plot(t_vec,GA_ti_32_slow(2,:),'Color',col_palette(4,:),'LineWidth',3)
% xlabel('time (s)')
% ylabel('coherence')
% plot(timevec_stats(statT32.mask),statT32.mask(statT32.mask).*0.04,'Color',col_palette(1,:),'LineWidth',2)
% plot(timevec_stats(statD32.mask),statD32.mask(statD32.mask).*0.0425,'Color',col_palette(4,:),'LineWidth',2)
% ylabel('coherence')
% xlabel('time (s)')
% xlim([-0.2 0.5])
% xticks([0:0.5:0.5])
% ylim([0.01 0.045])
% yticks([0:0.02:0.04])
% box off
% legend('guided T fast', 'guided T slow','guided D fast', 'guided D slow',Location='southeast')
% 
% 
% print(fig,fullfile(cohfigpth,'coh_fast_slow'),'-dpng','-r0')
% print(fig,fullfile(cohfigpth,'coh_fast_slow'),'-dsvg','-r0')
% 
% 
% fig = figure('Position',[0 0  1980 1080/2.5]);
% 
% timevec_stats = timevec(500:end);
% 
% 
% %% Unguided
% % for set size 32
% u16f = cell(1,length(subj));
% u16s = cell(1,length(subj));
% u32f = cell(1,length(subj));
% u32s = cell(1,length(subj));
% 
% for s = 1:length(subj)
% 
%     % 16
% 
%     % fast
%     u = squeeze(((niT(1,s,1,:) + niD(1,s,1,:)))./2);
% 
%     ERF.avg(1,:) = u;
%     u16f{s} = ERF;
% 
%     % slow
%     u = squeeze(((niT(1,s,2,:) + niD(1,s,2,:)))./2);
%     ERF.avg(1,:) = u;
%     u16s{s} = ERF;
% 
% 
%     % 32
% 
%     % fast
%     u = squeeze(((niT(2,s,1,:) + niD(2,s,1,:)))./2);
% 
%     ERF.avg(1,:) = u;
%     u32f{s} = ERF;
% 
%     % slow
%     u = squeeze(((niT(2,s,2,:) + niD(2,s,2,:)))./2);
%     ERF.avg(1,:) = D;
%     u32s{s} = ERF;
% 
% 
% end
% 
% cfg         = [];
% 
% cfg.method           = 'montecarlo';                                            % cluster-based
% cfg.channel          = ERF.label{1};
% cfg.statistic        = 'depsamplesT';                                           % within subject
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
% cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution   
% cfg.minnbchan        = 0;
% cfg.alpha            = 0.05;
% cfg.numrandomization = 500;
% cfg.latency = [0 min_rt];
% Nsubj  = length(subj);
% design = zeros(2, Nsubj*2);
% design(1,:) = [1:Nsubj 1:Nsubj];
% design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
% 
% cfg.design = design;
% cfg.uvar   = 1;
% cfg.ivar   = 2;
% 
% % compare Target and Distractor fast vs slow
% 
% cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
% cfg.clustertail = 1;
% 
% [statu16] = ft_timelockstatistics(cfg, u16f{:}, u16s{:});
% 
% cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
% cfg.clustertail = -1;
% [statu32] = ft_timelockstatistics(cfg, u32f{:}, u32s{:});
% 
% fig = figure('Position',[0 0  1980 1080/2.5]);
% 
% % set size 16
% ax1 = subplot(131);
% plot(t_vec,mean(GA_ni_16_fast,1),'Color',col_palette(5,:),'LineWidth',3)
% hold on
% plot(t_vec,mean(GA_ni_16_slow,1),'Color',col_palette(6,:),'LineWidth',3)
% 
% ylabel('coherence')
% xlabel('time (s)')
% xlim([-0.1 0.5])
% xticks([0:0.5:0.5])
% ylim([0.01 0.05])
% yticks([0:0.025:0.05])
% box off
% 
% legend('guided D fast', 'guided D slow',Location='southeast')
% 
% 
% % set size 32
% ax2 = subplot(132);
% plot(t_vec,mean(GA_ni_32_fast,1),'Color',col_palette(5,:),'LineWidth',3)
% hold on
% plot(t_vec,mean(GA_ni_32_slow,1),'Color',col_palette(6,:),'LineWidth',3)
% 
% 
% ylabel('coherence')
% xlabel('time (s)')
% xlim([-0.1 0.5])
% xticks([0:0.5:0.5])
% ylim([0.01 0.05])
% yticks([0:0.025:0.05])
% box off
% legend('guided D fast', 'guided D slow',Location='southeast')
% 
% 
% print(fig,fullfile(cohfigpth,'coh_fast_slow_ung'),'-dpng','-r0')
% print(fig,fullfile(cohfigpth,'coh_fast_slow_ung'),'-dsvg','-r0')
% 
