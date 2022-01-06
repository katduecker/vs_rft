%% VS + RFT
% PhD project 2

% coherence using hilbert transform on maxfiltered data

% [c] Katharina Duecker

% 1. Load in data
% 2. Split up based on Target tagging frequency

%clear all; close all; clc; beep off;

function a_coh_hilb(s)
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');

% matlab scripts path
mtpth = fullfile(pth,'matlab scripts','tfrs');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb');
cohsoipth = fullfile(cohpth,'soi');
cohrespth = fullfile(cohpth,'coh');

mkdir(cohrespth)

ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% load documentation: subjects with identified soi
load(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_soi_coh.mat'))

keep_subj = ~strcmp(mergesubj(:,5), 'no soi');

subjfolds = subjfolds(keep_subj);

% load soi
load(fullfile(cohsoipth,subjfolds{s},'soi_stat.mat'))
senstype = str2double(cellfun(@(x) x(end),soi_stat,'UniformOutput',false));
soimag = soi_stat(senstype == 1);
soigrad = soi_stat(senstype~=1);
save(fullfile(cohsoipth,subjfolds{s},'soi_stat.mat'),'soi_stat','soigrad','soimag')

% load trial structure per subject

load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'))

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,2)-fs,meginfo.alltrl_bl{p}(:,3)+fs*2,repmat(-2.5*fs,size(meginfo.alltrl_bl{p},1),1)];
    trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;
end


% read in artefactual sensors
fid = fopen(fullfile(dtpth,subjfolds{s},'artef_sens.txt'),'r');
sc = 1;
noisy_sens = {};
sens = fgetl(fid);
while isstr(sens)
    noisy_sens{sc} = sens; 
    sc = sc+1;
    sens = fgetl(fid);
end
fclose(fid)

% list fif files
d = dir(fullfile(dtpth,subjfolds{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load in maxfiltered data using trial structure
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = {'MEG'};
    
    % load data and diode separately (to be able to reject ICA)
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);
    % diodes
    cfg.channel = {'MISC004', 'MISC005'};
    diodes_trl{p}= ft_preprocessing(cfg);
end

data = ft_appenddata([],dtprt{:});
diodes = ft_appenddata([], diodes_trl{:});

% fix sampleinfo
sinfo = zeros(length(data.trial),1);
% first trial samples: 0 to length of data in samples
sinfo(1,2) = size(data.trial{1},2);
for t = 2:length(data.trial)
    % next trial: one sample after end of previous
    sinfo(t,1) = sinfo(t-1,2) + 1;
    % start sample + length of trial - 1
    sinfo(t,2) = sinfo(t,1) + size(data.trial{t},2);
end

data.sampleinfo = sinfo;

% select trials to be kept
if ~isempty(meginfo.rejtrl_meg)
    cfg = [];
    cfg.trials = ~ismember(1:length(data.sampleinfo),meginfo.rejtrl_meg);
    data = ft_selectdata(cfg,data);  
    diodes = ft_selectdata(cfg,diodes);

    meginfo.alltrl_list(meginfo.rejtrl_meg,:) = [];
end

if ~isempty(idx_rt_diff) && ~ignore_diff
    cfg = [];
    cfg.trials = ~ismember(1:length(data.sampleinfo),idx_rt_diff);
    data = ft_selectdata(cfg,data); 
    diodes = ft_selectdata(cfg,diodes);
    meginfo.alltrl_list(idx_rt_diff,:) = [];

end


% load ICA data and components
load(fullfile(icapth, [subjfolds{s},'_ica.mat']))

cfg = [];
cfg.component = badcomps; % to be removed component(s)
dataclean = ft_rejectcomponent(cfg, dataICA, data);
clear dataICA

% average grad positions
d = dir(fullfile(megpth,subjfolds{s}, 'meg'));
files = {d.name};
files = files(strncmp(files,'part',4));

% average channel positions
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

cfg = [];
cfg.resamplefs = 250;
dataresamp = ft_resampledata(cfg,data);

%% Coherence for tagging frequency

load(fullfile(pth,'experiment','trigdef.mat'))

% select trials 60 Hz Target

% find triggers
% T color: yellow (misc 4); T freq: 60 Hz
trig60Ty = [trigdef{cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'126067'),'UniformOutput',false)),1}];
% T color: cyan (misc 5); T freq: 60 Hz
trig60Tc = [trigdef{cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'216067'),'UniformOutput',false)),1}];
% T color: yellow; T freq: 67 Hz
trig67Ty = [trigdef{cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'126760'),'UniformOutput',false)),1}];
% T color: cyan; T freq: 67 Hz
trig67Tc = [trigdef{cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'216760'),'UniformOutput',false)),1}];

% find trials
trl60Ty = ismember([meginfo.alltrl_list(:,1)],trig60Ty);
trl60Tc = ismember([meginfo.alltrl_list(:,1)],trig60Tc);
trl67Ty = ismember([meginfo.alltrl_list(:,1)],trig67Ty);
trl67Tc = ismember([meginfo.alltrl_list(:,1)],trig67Tc);

% select data
cfg = [];
cfg.trials = trl60Ty;
data60Tmisc4 = ft_selectdata(cfg, data);
cfg.trials = trl60Tc;
data60Tmisc5 = ft_selectdata(cfg, data);
cfg.trials = trl67Ty;
data67Tmisc4 = ft_selectdata(cfg, data);
cfg.trials = trl67Tc;
data67Tmisc5 = ft_selectdata(cfg, data);

% sanity
length(data.trial) == (length(data60Tmisc4.trial) + length(data67Tmisc4.trial) + length(data60Tmisc5.trial) + length(data67Tmisc5.trial))

% replace diode labels with "diode T" (target freq) and "diode D"
% (disctractor freq)

% select T
cfg = [];
cfg.channel = 'MISC004';
data60Tmisc4_T = ft_selectdata(cfg,data60Tmisc4);
data60Tmisc4_T.label = {'diode T'};
data67Tmisc4_T = ft_selectdata(cfg,data67Tmisc4);
data67Tmisc4_T.label = {'diode T'};
cfg.channel = 'MISC005';
data60Tmisc5_T = ft_selectdata(cfg,data60Tmisc5);
data60Tmisc5_T.label = {'diode T'};
data67Tmisc5_T = ft_selectdata(cfg,data67Tmisc5);
data67Tmisc5_T.label = {'diode T'};

% select D
cfg = [];
cfg.channel = 'MISC005';
data60Tmisc4_D = ft_selectdata(cfg,data60Tmisc4);
data60Tmisc4_D.label = {'diode D'};
data67Tmisc4_D = ft_selectdata(cfg,data67Tmisc4);
data67Tmisc4_D.label = {'diode D'};
cfg.channel = 'MISC004';
data60Tmisc5_D = ft_selectdata(cfg,data60Tmisc5);
data60Tmisc5_D.label = {'diode D'};
data67Tmisc5_D = ft_selectdata(cfg,data67Tmisc5);
data67Tmisc5_D.label = {'diode D'};

% select only MEG channels
cfg = [];
cfg.channel = soi_stat;
data60Tmisc4_M = ft_selectdata(cfg,data60Tmisc4);
data60Tmisc5_M = ft_selectdata(cfg,data60Tmisc5);
data67Tmisc4_M = ft_selectdata(cfg,data67Tmisc4);
data67Tmisc5_M = ft_selectdata(cfg,data67Tmisc5);

% append with target and distractor label
data60T4 = ft_appenddata([],data60Tmisc4_M,data60Tmisc4_T,data60Tmisc4_D);
data60T5 = ft_appenddata([],data60Tmisc5_M,data60Tmisc5_T,data60Tmisc5_D);
data67T4 = ft_appenddata([],data67Tmisc4_M,data67Tmisc4_T,data67Tmisc4_D);
data67T5 = ft_appenddata([],data67Tmisc5_M,data67Tmisc5_T,data67Tmisc5_D);


% append data sets based on target freq
data60T = ft_appenddata([],data60T4,data60T5);
data67T = ft_appenddata([],data67T4,data67T5);

clear data60Tm* data67Tm* data

clear data60T4 data60T5 data67T4 data67T5

data60T.fsample = 1000;
data67T.fsample = 1000;

%% Coherence
foi = 55:75;
frqwdth = 2;

% Target freq 60 Hz
% gradiometers
[coh60.coh60Tgrad, psd60.psdmeg60Tgrad, psd60.psdmisc60Tgrad, csd60.csd60Tgrad] = kd_coh_hilb_fun(data60T,'diode T', 'MEGGRAD',foi, frqwdth);
% coherence distractor
[coh60.coh60Dgrad, psd60.psdmeg60Dgrad, psd60.psdmisc60Dgrad, csd60.csd60Dgrad] = kd_coh_hilb_fun(data60T,'diode D', 'MEGGRAD', foi, frqwdth);

% magnetometers
[coh60.coh60Tmag, psd60.psdmeg60Tmag, psd60.psdmisc60Tmag, csd60.csd60Tmag] = kd_coh_hilb_fun(data60T,'diode T', 'MEGMAG',foi, frqwdth);
% coherence distractor
[coh60.coh60Dmag, psd60.psdmeg60Dmag, psd60.psdmisc60Dmag, csd60.csd60Dmag] = kd_coh_hilb_fun(data60T,'diode D', 'MEGMAG', foi, frqwdth);


mkdir(fullfile(cohrespth,subjfolds{s}))
save(fullfile(cohrespth,subjfolds{s},['coh_freqw_',num2str(frqwdth),'_T60.mat']),'coh60','psd60','csd60')

clear coh60 psd60 csd60

% Target freq 67 Hz
% gradiometers
[coh67.coh67Tgrad, psd67.psdmeg67Tgrad, psd67.psdmisc67Tgrad, csd67.csd67Tgrad] = kd_coh_hilb_fun(data67T,'diode T', 'MEGGRAD',foi, frqwdth);
% coherence distractor
[coh67.coh67Dgrad, psd67.psdmeg67Dgrad, psd67.psdmisc67Dgrad, csd67.csd67Dgrad] = kd_coh_hilb_fun(data67T,'diode D', 'MEGGRAD', foi, frqwdth);

% magnetometers
[coh67.coh67Tmag, psd67.psdmeg67Tmag, psd67.psdmisc67Tmag, csd67.csd67Tmag] = kd_coh_hilb_fun(data67T,'diode T', 'MEGMAG',foi, frqwdth);
% coherence distractor
[coh67.coh67Dmag, psd67.psdmeg67Dmag, psd67.psdmisc67Dmag, csd67.csd67Dmag] = kd_coh_hilb_fun(data67T,'diode D', 'MEGMAG', foi, frqwdth);


save(fullfile(cohrespth,subjfolds{s},['coh_freqw_',num2str(frqwdth),'_T67.mat']),'coh67','psd67','csd67')
