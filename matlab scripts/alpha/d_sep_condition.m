%% VS + RFT
% PhD project 2

% d. contrast conditions

% Inputs
% - s: subject index
% - winl: window length in seconds
% - c_idx: condition index

% Output
% - TFR contrasts between conditions

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 28/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions

function d_sep_condition(s,winl,c_idx)

condi = {{'ni'},{'ti'},{'16t'},{'32t'}};

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha','pow');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi/');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

d = dir(cohpth);
d = {d.name};
folds = d(strncmp(d,'202',3));

load('alpha_align_vec.mat')
folds(excl_alpha) = [];

load(fullfile(mergepth, folds{s},'trl_overlap_meg_el_rsp.mat'))
load(fullfile(pth, 'experiment','trigdef.mat'))

trl = rspinfo.trl(rspinfo.keeptrl_rsp,:);
trl = trl(meginfo.keeptrl_all,:);

cur_cond = condi{c_idx};

% find triggers associated with current condition
for c = 1:length(cur_cond)
    cond_idx(c,:) = cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),cur_cond{c}),'UniformOutput',false));
end
cond_idx = sum(cond_idx,1) == length(cur_cond);

%% delete file if it exists

if exist(fullfile(outpth,folds{s},[strjoin(cur_cond,'_'),'_',num2str(winl*10),'.mat']),'file')
    delete(fullfile(outpth,folds{s},[strjoin(cur_cond,'_'),'_',num2str(winl*10),'.mat']))
end

% find trials associated with current condition
trig_oi = trigdef(cond_idx,1);

trl_idx = ismember([trl{:,1}],[trig_oi{:}]);
trl_oi = trl(ismember([trl{:,1}],[trig_oi{:}]),:);

load(fullfile(outpth,folds{s},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha')

cfg = [];
cfg.trials = find(trl_idx);
cfg.avgoverrpt = 'yes';
TFR_cond = ft_selectdata(cfg,TFR_alpha);

powspct = TFR_cond.powspctrm;

save(fullfile(outpth,folds{s},[strjoin(cur_cond,'_'),'_',num2str(winl*10),'.mat']),'powspct','trl_idx','trl_oi')