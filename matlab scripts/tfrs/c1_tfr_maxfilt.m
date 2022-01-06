%% VS + RFT
% PhD project 2

% TFR of power on maxfiltered data

% [c] Katharina Duecker

% 1. TFR on ICA whole data
% 2. Hanning 1 second

%clear all; close all; clc; beep off

function c1_tfr_maxfilt(s)

% window length
winl = 1;
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');
% save tfr
tfrpth = fullfile(pth,'results','meg', '4 TFR power');
% figure paths
figpth = fullfile(pth,'results','meg', 'figures', 'TFR power');
mkdir(tfrpth)
mkdir(figpth)
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% identify bad components
%for s = 1:length(subjfolds)

load('soi_tfr_subj.mat')
load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'alltrl_bl')

% trial structure to load in trl
for p = 1:length(alltrl_bl)
    trlstruct{p} = [alltrl_bl{p}(:,2)-fs,alltrl_bl{p}(:,3)+fs,repmat(-2.5*fs,size(alltrl_bl{p},1),1)];
end

% load in maxfiltered data using trial structure

% list fif files
d = dir(fullfile(dtpth,subjfolds{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = {'MEG'};
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);
    % diodes
    cfg.channel = {'MISC004', 'MISC005'};
    diodes_trl{p}= ft_preprocessing(cfg);
end

data = ft_appenddata([],dtprt{:});
diodes = ft_appenddata([], diodes_trl{:});

% load ICA data and components
load(fullfile(icapth, [subjfolds{s},'_ica.mat']))

cfg = [];
cfg.component = badcomps; % to be removed component(s)
dataclean = ft_rejectcomponent(cfg, dataICA, data);
clear dataICA

% % compare before and after rejection
% cfg = [];
% cfg.layout = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
% cfg.viewmode = 'vertical';
% ft_databrowser(cfg, data)
% ft_databrowser(cfg, dataclean)
% clear data

% average grad positions
d = dir(fullfile(megpth,subjfolds{s}, 'meg'));
files = {d.name};
files = files(strncmp(files,'part',4));

grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(megpth,subjfolds{s},'meg',files{fl}))];
end

% mean channel position
mGrad = grad(1);
% average grad structure
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
clear grad

% concatenate cleaned data and diodes
data = ft_appenddata([],dataclean,diodes);
clear diodes diodes_trl dataclean
%% TFR all trials

cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
% cfg.channel = {'MEG','MISC004','MISC005'};
cfg.pad = 10;
cfg.padtype = 'zero';
cfg.foi = 54:1:80;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
cfg.toi = -1.25:0.05:3.5;
cfg.keeptrials = 'no';
TFRall = ft_freqanalysis(cfg,data);
TFRall.grad = mGrad;

save(fullfile(tfrpth,[subjfolds{s},'_TFRalltrl.mat']),'TFRall')
