%% VS + RFT0)
% PhD project 2

% g. Split trials into conditions

% [c] Katharina Duecker

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions


function g_split_trl(s,rej_sac)
%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path
ldiode = {'MISC004','MISC005'}; % diode label
rft_freq = [60,67];             % tagging frequencies
phshft = [0,0.65];              % phase shift per tagging frequency
fs = 1000;                      % sampling rate

dtpth = fullfile(pth, 'data'); % data path
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))
ft_defaults;

dtpth = fullfile(pth,'data');                                       % raw data
maxfpth = fullfile(pth,'results','meg', '1 maxfilter');             % max filter
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');   % scripts
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
icapth = fullfile(pth,'results','meg', '3 ICA', '1 all subj');
respth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid');

%% load in data and discard artefactual trials
load(fullfile(scriptpth,'idx_subjoi.mat'))                          % usable subjects index

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

subj = subjfolds(usable_idx);

% load in trial structure
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*2.5,meginfo.alltrl_bl{p}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];
end

% list maxfiltered data
d = dir(fullfile(maxfpth,subj{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load data
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(maxfpth,subj{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = 'MEG';
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);
    % diodes
    cfg.channel = {'MISC004', 'MISC005'};
    diodes_trl{p}= ft_preprocessing(cfg);
end

data = ft_appenddata([],dtprt{:});
data.hdr = dtprt{1}.hdr;
diodes = ft_appenddata([], diodes_trl{:});

meginfo.alltrl_list = meginfo.alltrl_list(logical(meginfo.keeptrl_all),:);

% discard strange trials
cfg = [];
cfg.trials =find(meginfo.keeptrl_all);                   % select clean trials
data = ft_selectdata(cfg,data);
diodes = ft_selectdata(cfg,diodes);

% sampleinfo
trl_end_samp = cell2mat(cellfun(@length,data.trial,'UniformOutput',false)).*[1:length(data.trial)];
data.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];
diodes.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];

%% load in ICA and reject comps

% load ICA data and components
load(fullfile(icapth, [subj{s},'_ica.mat']))

cfg = [];
cfg.component = badcomps; % to be removed component(s)
dataclean = ft_rejectcomponent(cfg, dataICA, data);
clear dataICA

% average grad positions
d = dir(fullfile(dtpth,subj{s}, 'meg'));
files = {d.name};
files = files(strncmp(files,'part',4));

%% average channel positions
grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(dtpth,subj{s},'meg',files{fl}))];
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
data.grad = mGrad;

%% reject trials containing eye movement (?)
if rej_sac
    trl_keep = 1:length(data.trial);
    trl_keep = ~ismember(trl_keep,elinfo.sactrial_lib);
    meginfo.alltrl_list = meginfo.alltrl_list(trl_keep,:);

    cfg = [];
    cfg.trials = trl_keep;
    data = ft_selectdata(cfg,data);
    
end

%% Extract data with current trigger code

load(fullfile(pth,'experiment','trigdef.mat'))

for trg_idx = 5:length(trigdef)
    % find trials
    trlcur = ismember(meginfo.alltrl_list(:,1),trigdef{trg_idx});

    % select data
    cfg = [];
    cfg.trials = trlcur;
    data_cur = ft_selectdata(cfg, data);

    % select T: Target color: end-5: channel index (1: MISC004, 2: MISC005)
    cfg = [];
    chanT = ldiode{str2num(trigdef{trg_idx,2}(end-5))};
    cfg.channel = chanT;
    data_T = ft_selectdata(cfg,data_cur);
    data_T.label = {'diode T'};
    
    % replace with perfect sinusoid
    freqT = str2num(trigdef{trg_idx,2}(end-3:end-2));           % target frequency
    data_T = kd_replace_diode_sinu(data_T,1,freqT,phshft(rft_freq == freqT),-2.5,fs);
    % phase shift: 0 for 60 Hz, 0.65 for 67
    
    % select D (the channel that was not selected above)
    cfg = [];
    cfg.channel = ldiode{~strcmp(ldiode,chanT)};
    data_D = ft_selectdata(cfg,data_cur);
    data_D.label = {'diode D'};
    
    % replace with perfect sinusoid
    freqD = str2num(trigdef{trg_idx,2}(end-1:end));
    data_D = kd_replace_diode_sinu(data_D,1,freqD,phshft(rft_freq == freqD),-2.5,fs);
    
    % only select MEG
    cfg= [];
    cfg.channel = 'MEG';
    data_cur = ft_selectdata(cfg,data_cur);


    % append with target and distractor label
    data_trig = ft_appenddata([],data_cur,data_T,data_D);
    data_trig.fsample = 1000;

    mkdir(fullfile(respth,subj{s}))
    if rej_sac
        mkdir(fullfile(respth,subj{s}),'rej_sac')
        save(fullfile(respth,subj{s},'rej_sac',['data_',trigdef{trg_idx,2},'.mat']),'data_trig')
    else
        save(fullfile(respth,subj{s},['data_',trigdef{trg_idx,2},'.mat']),'data_trig')
    end
    clear data_D data_T data_cur data_trig
end
