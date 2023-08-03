%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d1. Coherence for alpha high/low

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023
 
% Input:
% -s: subject index
% - c_idx: condition index
% - fwdth: bandwith of bpfilter applied before Hilbert transform
% (recommended: 5 Hz)
% - toi: time interval for alpha median split
% - split_ta_tp: separate for target absent/present? if 0, the trials will
% still be split into ta/tp and then recombined

% Output
% Coherence for high low separated for each condition: coh6067_high/...low & coh6760_high/...low

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)



function d1_rift_alpha_split(s, c_idx,fwdth,toi,split_ta_tp)

% Input
% s: subject index
% c_idx: condition index (1-4)
% fwdth : width of bp filter before hilbert transform (in Hz)
% toi: time window used to split alpha power
% ta_tp: split into target absent/present?

%% Paths
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

addpath(fullfile(pth,'matlab scripts','alpha'))


inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');

alphapth = fullfile(pth,'results','meg','6 Alpha');
alphapowpth = fullfile(alphapth,'not align','pow align iaf');
soipth = fullfile(alphapth,'iaf_soi');
cohsoipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid/');

rtpth = fullfile(pth,'results','behavior','alpha');
mkdir(fullfile(rtpth))

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','alpha high low');

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

condi_all = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

condi = condi_all{c_idx};
% load(fullfile(alphapowpth,subj{s},'data_winl_5.mat'),'TFR_alpha','perf_TFR')
% load SOI
load(fullfile(soipth,subj{s},'iaf_soi.mat'))

% split into Target absent/Target present?

if split_ta_tp
    % if yes, separate trials for ta/tp
    ta_tp = {'ta','tp'};
else
    % if not, just concatenate a t to the set size (in filenames) which
    % looks for ta/tp
    ta_tp = {'t'};
end

d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];


for ti = 1:length(ta_tp)

 % find relevant data files
condi_files = zeros(length(files),1);
for c = 1:length(condi)
    
    condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c}),'UniformOutput',false))';
    
end

condi_files = condi_files == length(condi);   
condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,[condi{2}(end-2:end-1),ta_tp{ti}]),'UniformOutput',false))';


files6067 =  condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'6067'),'UniformOutput',false))';
files6760 =  condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'6760'),'UniformOutput',false))';

files6067 = files6067 == 3;
files6760 = files6760 == 3;


% load files 6067
files6067 = files(files6067);

load(fullfile(inpth,subj{s},files6067{1}));
trl_idx = find(trlcur);

% performance in first trials
trl_perf6067 = perf_cur;

% load data
data_load = data_trig;

% data_load = rmfield(data_load,'sampleinfo');
% load & append data
for f = 2:length(files6067)
    clear data_trig trlcur perf_cur
    load(fullfile(inpth,subj{s},files6067{f}));
    % append data
    data_load = ft_appenddata([],data_load,data_trig);
    %append performance
    trl_perf6067 = [trl_perf6067;perf_cur];
end

data6067 = data_load;

clear data_load


% load files 6760
files6760 = files(files6760);

load(fullfile(inpth,subj{s},files6760{1}));
trl_idx = find(trlcur);

% performance in first trials
trl_perf6760 = perf_cur;

% load data
data_load = data_trig;
% data_load = rmfield(data_load,'sampleinfo');
% load & append data
for f = 2:length(files6760)
    clear data_trig trlcur perf_cur
    load(fullfile(inpth,subj{s},files6760{f}));
    % append data
    data_load = ft_appenddata([],data_load,data_trig);
    %append performance
    trl_perf6760 = [trl_perf6760;perf_cur];
end

data6760 = data_load;

clear data_load



% trim based on reaction time
all_rt = [trl_perf6067{:,3},trl_perf6760{:,3}];
m_rt = mean(all_rt);
std_rt = std(all_rt);


trl_rej6067 = ([trl_perf6067{:,3}] < m_rt - std_rt*3)+([trl_perf6067{:,3}] > m_rt + std_rt*3);
trl_rej6760 = ([trl_perf6760{:,3}] < m_rt - std_rt*3)+([trl_perf6760{:,3}] > m_rt + std_rt*3);

% reject trials with outliers
cfg = [];
cfg.trials = ~trl_rej6067;
data6067 = ft_selectdata(cfg,data6067);
cfg.trials = ~trl_rej6760;
data6760 = ft_selectdata(cfg,data6760);

% TFR for condition (ta/tp)
data_condi = ft_appenddata([],data6067,data6760);


trl_perf6067 = trl_perf6067(~trl_rej6067,:);
trl_perf6760 = trl_perf6760(~trl_rej6760,:);


%% median split alpha power
winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = soi_grad;
cfg.taper = 'hanning';
cfg.foi = iaf_grad;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'yes';
TFR_alpha = ft_freqanalysis(cfg,data_condi);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);

% average over baseline and extract IAF power
cfg = [];

cfg.latency = toi;
cfg.avgovertime = 'yes';

cfg.channel = soi_grad_cmb;

cfg.avgoverchan = 'yes';
IAFpow = ft_selectdata(cfg,TFR_alpha);

iaf_pow = IAFpow.powspctrm;

m_iaf = median(iaf_pow);


clear data_condi TFR_alpha IAFpow

%% Performance for high and low

% select trials wrt tagging freq
% re-do TFR
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = soi_grad;
cfg.taper = 'hanning';
cfg.foi = iaf_grad;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'yes';
TFR6067 = ft_freqanalysis(cfg,data6067);
TFR6760 = ft_freqanalysis(cfg,data6760);

cfg = [];
cfg.method = 'sum';
TFR6067 = ft_combineplanar(cfg,TFR6067);
TFR6760 = ft_combineplanar(cfg,TFR6760);

% average over baseline and extract IAF power
cfg = [];

cfg.latency = toi;
cfg.avgovertime = 'yes';

cfg.channel = soi_grad_cmb;
cfg.avgoverchan = 'yes';
IAF6067 = ft_selectdata(cfg,TFR6067);
IAF6760 = ft_selectdata(cfg,TFR6760);

% trial indices
trl6067_high = find(IAF6067.powspctrm>m_iaf);
trl6067_low = find(IAF6067.powspctrm<m_iaf);

trl6760_high = find(IAF6760.powspctrm>m_iaf);
trl6760_low = find(IAF6760.powspctrm<m_iaf);


%% select performance
rt_high = mean([[trl_perf6067{trl6067_high,3}],[trl_perf6760{trl6760_high,3}]]);
rt_low = mean([[trl_perf6067{trl6067_low,3}],[trl_perf6760{trl6760_low,3}]]);


hit_high = (sum(strcmp(trl_perf6067(trl6067_high,2),'h')) + sum(strcmp(trl_perf6760(trl6760_high,2),'h')))/...
    (length(trl6067_high)+length(trl6760_high));

hit_low = (sum(strcmp(trl_perf6067(trl6067_low,2),'h')) + sum(strcmp(trl_perf6760(trl6760_low,2),'h')))/...
    (length(trl6067_low)+length(trl6760_low));


clear IAF6* TFR6*

%% Coherence

% select trials
cfg = [];
cfg.trials = trl6067_high;
data6067_high = ft_selectdata(cfg,data6067);

cfg.trials = trl6067_low;
data6067_low = ft_selectdata(cfg,data6067);

cfg.trials = trl6760_high;
data6760_high = ft_selectdata(cfg,data6760);

cfg.trials = trl6760_low;
data6760_low = ft_selectdata(cfg,data6760);


[coh6067_high.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_high,'diode T', {'MEGGRAD'},60, fwdth,{'but','twopass'});
[coh6067_high.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_high,'diode D', {'MEGGRAD'},67, fwdth,{'but','twopass'});
[coh6760_high.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_high,'diode T', {'MEGGRAD'},67, fwdth,{'but','twopass'});
[coh6760_high.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_high,'diode D', {'MEGGRAD'},60, fwdth,{'but','twopass'});


[coh6067_low.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_low,'diode T', {'MEGGRAD'},60, fwdth,{'but','twopass'});
[coh6067_low.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_low,'diode D', {'MEGGRAD'},67, fwdth,{'but','twopass'});
[coh6760_low.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_low,'diode T', {'MEGGRAD'},67, fwdth,{'but','twopass'});
[coh6760_low.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_low,'diode D', {'MEGGRAD'},60, fwdth,{'but','twopass'});

% baseline correct

% high

% 6067
bsl = mean(coh6067_high.cohTgrad(:,:,1000:2500),3);
coh6067_high.bslcor.cohTgrad = coh6067_high.cohTgrad - bsl;

bsl = mean(coh6067_high.cohDgrad(:,:,1000:2500),3);
coh6067_high.bslcor.cohDgrad = coh6067_high.cohDgrad - bsl;

% 6760
bsl = mean(coh6760_high.cohTgrad(:,:,1000:2500),3);
coh6760_high.bslcor.cohTgrad = coh6760_high.cohTgrad - bsl;

bsl = mean(coh6760_high.cohDgrad(:,:,1000:2500),3);
coh6760_high.bslcor.cohDgrad = coh6760_high.cohDgrad - bsl;

% low

% 6067
bsl = mean(coh6067_low.cohTgrad(:,:,1000:2500),3);
coh6067_low.bslcor.cohTgrad = coh6067_low.cohTgrad - bsl;

bsl = mean(coh6067_low.cohDgrad(:,:,1000:2500),3);
coh6067_low.bslcor.cohDgrad = coh6067_low.cohDgrad - bsl;

% 6760
bsl = mean(coh6760_low.cohTgrad(:,:,1000:2500),3);
coh6760_low.bslcor.cohTgrad = coh6760_low.cohTgrad - bsl;

bsl = mean(coh6760_low.cohDgrad(:,:,1000:2500),3);
coh6760_low.bslcor.cohDgrad = coh6760_low.cohDgrad - bsl;

% save coherence results
mkdir(fullfile(cohpth,subj{s}))


toi_str = arrayfun(@num2str,toi.*1000,'UniformOutput',false);

condname = strjoin(condi,'_');
save(fullfile(cohpth,subj{s},[condname(1:end-1),ta_tp{ti},'_fwidth',num2str(fwdth),strjoin(toi_str,'_'),'_condi.mat']),'coh6067_high','coh6067_low','coh6760_high','coh6760_low','trl6067_high','trl6067_low','trl6760_high','trl6760_low','iaf_pow','rt_high','rt_low','hit_high','hit_low')
end

