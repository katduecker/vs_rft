%% VS + RFT
% PhD project 2

% coherence using hilbert transform on maxfiltered data
% run coherence for defined condition

% [c] Katharina Duecker

% 1. Load in data
% 2. Split up based on Target tagging frequency

%clear all; close all; clc; beep off;

function d_extract_data_coh(s,trigcode)

% s: subject ID
% trigcode: trigger code as defined in trigdef

% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');

ldiode = {'MISC004','MISC005'};
% matlab scripts path
mtpth = fullfile(pth,'matlab scripts','tfrs');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb');
cohsoipth = fullfile(cohpth,'soi');
cohdatpth = fullfile(cohpth,'data sep');

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

% make output directory
mkdir(fullfile(cohdatpth,subjfolds{s}))

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
cfg.resamplefs = 500;
dataresamp = ft_resampledata(cfg,data);

%% Coherence for tagging frequency

load(fullfile(pth,'experiment','trigdef.mat'))

% select trials 60 Hz Target

% % find triggers
% for c = 1:length(colfreqcode)
%     % separation in frequency and T color always needed
%     trigcolf(c,:) = [trigdef{cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),colfreqcode{c}),'UniformOutput',false)),1}];
%     % if an additional group name is defined, separate into that as well
%     if ~isempty(groupname)
%          triggroup = [trigdef{cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),groupname),'UniformOutput',false)),1}];
%          trigcombi(c,:) = intersect(trigcolf(c,:),triggroup);
%     else 
%         trigcombi = trigcolf;
%     end
% end

% find trials
trlcur = ismember([meginfo.alltrl_list(:,1)],trigcode);

% select data
cfg = [];
cfg.trials = trlcur;
data_cur = ft_selectdata(cfg, data);

% select T: Target color: 
cfg = [];
chanT = ldiode{str2num(trigdef{trigcode,2}(end-5))};
cfg.channel = chanT;
data_T = ft_selectdata(cfg,data_cur);
data_T.label = {'diode T'};

% select D (the channel that was not selected above)
cfg = [];
cfg.channel = ldiode{~strcmp(ldiode,chanT)};
data_D = ft_selectdata(cfg,data_cur);
data_D.label = {'diode D'};


% append with target and distractor label
data = ft_appenddata([],data_cur,data_T,data_D);
data.fsample = 1000;

save(fullfile(cohdatpth,subjfolds{s},['data_',trigdef{trigcode,2},'.mat']),'data')
