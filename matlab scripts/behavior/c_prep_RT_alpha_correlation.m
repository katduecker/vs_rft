%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c. prep confirmatory correlation: extract alpha power and RT for each
% trial

% Input
% -s : subject index


% Output
% csv file per participant, containing RT and alpha power for each trial 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Behavioural analyses
% a: performance per condition
% b. performance for alpha high vs low
% c: confirmatory analysis for median split sceptics: correlation alpha
% power ~ RT over trials


function c_prep_RT_correlation(s)
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','no noise diode');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha');
soipth = fullfile(alphapth,'iaf_soi');
outpth = fullfile(pth,'results','behavior','alpha','rt_alpha_regr');
mkdir(fullfile(outpth,'csv'))
mkdir(fullfile(outpth,'mat'))

%% load data condition
d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];

data_load = cell(1,length(files));
trl_load = {};
cfg = [];
for f = 1:length(files)
    
    load(fullfile(inpth,subj{s},files{f}))
    data_load{f} = data_trig;
    trl_load = [trl_load;perf_cur];
    
    clear data_trig perf_cur
end

data = ft_appenddata([],data_load{:});

clear data_load


% extract specs condition
load(fullfile(pth,'experiment','trigdef.mat'))

condi_specs = {'ti','32t','tp'};

rt_perf_alpha_condi = zeros(length(trl_load),4);

% store rt
rt_perf_alpha_condi(:,1) = [trl_load{:,3}];

% find trigger condition
% store: 1: guided; 1: set 32; 1: target present

for c = 1:length(condi_specs)
    trig = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x), strfind(trigdef(:,2),condi_specs{c}),'UniformOutput',false)),1));
       
    idx = ismember([trl_load{:,1}],trig);
    rt_perf_alpha_condi(:,c+2) = idx;
end

condi_comb = [[0,0,0];[0,0,1];[0,1,0];[0,1,1];[1,0,0];[1,0,1];[1,1,0];[1,1,1]];


h = 1;
condi_trl = zeros(length(condi_comb),2);

rej_trl = [];
for c = 1:length(condi_comb)
    h_start = h;
    while  h <= length(rt_perf_alpha_condi) && sum(rt_perf_alpha_condi(h,3:end) == condi_comb(c,:))==3
        h = h + 1;
    end
    h_end = h-1;
    condi_trl(c,:) = [h_start, h_end];
    
    m = mean([trl_load{condi_trl(1,1):condi_trl(1,2),3}]);
    std_cur = std([trl_load{condi_trl(1,1):condi_trl(1,2),3}]);
    
    rej_trl = [rej_trl,find(([trl_load{condi_trl(c,1):condi_trl(c,2),3}] < m-2*std_cur) + ([trl_load{condi_trl(c,1):condi_trl(c,2),3}] > m+2*std_cur)) + h_start-1];
end

keep_trl = ones(1,length(rt_perf_alpha_condi));
keep_trl(rej_trl) =0;

rt_perf_alpha_condi(rej_trl,:) = [];

% now add alpha power
load(fullfile(soipth,subj{s},'iaf_soi_not_align.mat'))

% TFR
winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = soi_grad;
cfg.taper = 'hanning';
cfg.foi = iaf_grad;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'yes';
cfg.trials = find(keep_trl);
TFR_alpha = ft_freqanalysis(cfg,data);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);

% average over time
cfg = [];
cfg.avgovertime = 'yes';
cfg.latency = [-1 0];

cfg.channel = soi_grad_cmb;
cfg.avgoverchan = 'yes';

IAF = ft_selectdata(cfg,TFR_alpha);

iaf_pow = IAF.powspctrm;
rt_log = (rt_perf_alpha_condi(:,1)-mean(rt_perf_alpha_condi(:,1)))./std(rt_perf_alpha_condi(:,1));
rt_perf_alpha_condi(:,2) = (iaf_pow - mean(iaf_pow))./std(iaf_pow);
rt_perf_alpha_condi(:,1) = rt_log;
T = array2table(rt_perf_alpha_condi,'VariableNames',{'RT','iafpow','gui','set32','target p'});



writetable(T,fullfile(outpth,'csv',[subj{s},'_rt_alpha.csv']))
save(fullfile(outpth,'mat',[subj{s},'_rt_alpha.mat']),'rt_perf_alpha_condi')