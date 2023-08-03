%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a. TFR of alpha over all trials

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs
% - s: subject index
% - c_idx: condition index (it's just 1 because we are looking at all data)
% - toi: start and end time point the sliding window is centered over
% - winl: window length
% Output
% - TFR of power with window length 1 sec, Hanning window, moved in 50 sec
% steps

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 23/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

function a_fun_alpha_pow(s,c_idx,toi,winl)

%condi = {{'data'},{'ni', '16t'},{'ni', '32t'},{'ti', '16t'},{'ti', '32t'}};

condi = {{'data'}};

condi = condi{c_idx};
%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha', 'pow');

rmpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

if exist(fullfile(outpth,[strjoin(condi,'_'),'_winl_',num2str(winl*10),'.mat']),'file')

    error('Alpha TFR exists for this subject.')
    
end

outpth = fullfile(outpth,subj{s});
mkdir(outpth)

% list condition files that contain 'condi'
d = dir(fullfile(inpth,subj{s}));

d = {d.name};
for c = 1:length(condi)
    cond_idx(c,:) = cell2mat(cellfun(@(x) ~isempty(x),regexp(d,condi(c)),'UniformOutput',false));
end

cond_idx = sum(cond_idx,1) == length(condi);
files = d(cond_idx);

load(fullfile(inpth,subj{s},files{1}));cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEG';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = toi(1):0.05:toi(2);
cfg.keeptrials = 'yes';
TFR_alpha = ft_freqanalysis(cfg,data);
data_load = data_trig;
perf_TFR = perf_cur;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,subj{s},files{f}));
    data_load = ft_appenddata([],data_load,data_trig);
    perf_TFR = [perf_TFR;perf_cur];
end

data = data_load;

data = rmfield(data,'cfg');
clear data_load

%% TFR
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEG';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = toi(1):0.05:toi(2);
cfg.keeptrials = 'yes';
TFR_alpha = ft_freqanalysis(cfg,data);

cfg = [];
cfg.avgoverrpt = 'yes';
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha);
save(fullfile(outpth,[strjoin(condi,'_'),'_winl_',num2str(winl*10),'.mat']),'TFR_alpha','TFR_alpha_avg','perf_TFR','-v7.3')